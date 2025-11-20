#!/bin/bash

aws_instance=$(aws ec2 describe-instances \
    --filters "Name=tag:Name,Values=Nextflow-runner" "Name=instance-state-name,Values=running" \
    --query "Reservations[0].Instances[0].InstanceId" \
    --output text)

if [ $? -ne 0 ] || [ -z "$aws_instance" ]; then
  echo "Error finding Nextflow-runner instance. Check your AWS configuration?"
  exit 1
fi

aws ssm start-session \
  --target "$aws_instance" \
  --document-name AWS-StartInteractiveCommand \
  --parameters '{"command":["sudo su - ec2-user"]}'
