#!/bin/bash
set -u

# Run the ScPCA-nf pipeline with options to specify the run mode and output
#
# Available RUN_MODE values are:
#   example:   run the example workflow
#   testing:      run the workflow in the development testing environment
#   staging:   run the full workflow with full data in the staging environment

GITHUB_TAG=${GITHUB_TAG:-main}
RUN_MODE=${RUN_MODE:-example}
RESUME=${RESUME:-false}
NEXTFLOW_PARAMS=${NEXTFLOW_PARAMS:-""}

date=$(date "+%Y-%m-%d")
datetime=$(date "+%Y-%m-%dT%H%M")
log_path=s3://ccdl-scpca-workdir-997241705947-us-east-1/logs/${RUN_MODE}/${date}

# Make sure environment includes local bin (where Nextflow is installed)
if ! [[ "$PATH" =~ "$HOME/.local/bin:$HOME/bin" ]]
then
    PATH="$HOME/.local/bin:$HOME/bin:$PATH"
fi
export PATH

# Get secrets from AWS Secrets Manager/1Password
AWS_SECRETS=$(aws secretsmanager get-secret-value --secret-id scpca_service_account_token | jq -r '.SecretString')
# AWS secrets are a key-value store: retrieve individual values with jq
OP_SERVICE_ACCOUNT_TOKEN=$(jq -r '.op_token' <<< $AWS_SECRETS)
export OP_SERVICE_ACCOUNT_TOKEN
TOWER_ACCESS_TOKEN=$(op read "$(jq -r '.op_seqera_token' <<< $AWS_SECRETS)")
export TOWER_ACCESS_TOKEN
TOWER_WORKSPACE_ID=$(op read "$(jq -r '.op_seqera_workspace' <<< $AWS_SECRETS)") # Use the OpenScPCA workspace
export TOWER_WORKSPACE_ID

# SLACK_WEBHOOK=$(op read "$(jq -r '.op_slack_webhook' <<< $AWS_SECRETS)")
# export SLACK_WEBHOOK

# slack_error() {
#   # function to create a slack message from an error log
#   log_file=$1
#   # add header and bullet points to the log file
#   message=$(printf "⚠️ Errors running scpca-nf pipeline:\n\n"; sed -e 's/^/• /' < "$log_file")
#   jq -n --arg message "$message" \
#     '{text: "Error running scpca-nf workflow.",
#       blocks: [{
#         type: "section",
#         text: {
#           type: "mrkdwn",
#           text: $message
#         }
#       }]
#     }' \
#     | curl --json @- "$SLACK_WEBHOOK"
# }

# move to nextflow app directory
cd /opt/nextflow
# create an empty log file to capture any errors
cat /dev/null > run_errors.log

# Define Nextflow profiles based on output mode
if [[ "$RUN_MODE" == "example" ]]; then
  profile="example"
elif [[ "$RUN_MODE" == "testing" ]]; then
  profile="ccdl_testing"
elif [[ "$RUN_MODE" == "staging" ]]; then
  profile="ccdl_staging"
else
  echo "Invalid RUN_MODE: $RUN_MODE. Must be one of: example, testing, staging." >> run_errors.log
  exit 1
fi

# add batch to profile
profile="${profile},batch"

if [ "$RESUME" == "true" ]; then
  resume_flag="-resume"
else
  resume_flag=""
fi

nextflow pull AlexsLemonade/scpca-nf -revision $GITHUB_TAG \
|| echo "Error pulling scpca-nf workflow with revision '$GITHUB_TAG'" >> run_errors.log

# post any errors from the workflow pull and exit
if [ -s run_errors.log ]; then
  # slack_error run_errors.log
  aws s3 cp run_errors.log "${log_path}/${datetime}_run_errors.log"
  exit 1
fi


nextflow run AlexsLemonade/scpca-nf \
  -revision $GITHUB_TAG \
  -profile $profile \
  -with-report "${datetime}_scpca_report.html" \
  -with-trace  "${datetime}_scpca_trace.txt" \
  -with-tower \
  $resume_flag \
  $NEXTFLOW_PARAMS \
  || echo "Error with scpca data run" >> run_errors.log

mv .nextflow.log "${datetime}_scpca-nf.log"


# replace any instances of TOWER_ACCESS_TOKEN in logs with masked value
sed -i "s/${TOWER_ACCESS_TOKEN}/<TOWER_ACCESS_TOKEN>/g" ./*.log*

# Copy logs to S3 and delete if successful
aws s3 cp . "${log_path}" \
  --recursive \
  --exclude "*" \
  --include "${datetime}_*" \
  && rm "${datetime}"_* \
  || echo "Error copying logs to S3" >> run_errors.log


# For example runs, compress the output directory and upload it to S3
# if [[ "$RUN_MODE" == "example" ]]; then
#   aws s3 cp --recursive "s3://scpca-nf-references/example-data/scpca_out" scpca_out \
#   && zip -r scpca_out.zip scpca_out \
#   && aws s3 cp scpca_out.zip "s3://scpca-nf-references/example-data/scpca_out.zip" \
#   || echo "Error uploading scpca_out.zip to S3" >> run_errors.log
#   rm -rf scpca_out scpca_out.zip
# fi

# Post any errors to slack
if [ -s run_errors.log ]; then
  # slack_error run_errors.log
  exit 1
else
  exit 0
fi
