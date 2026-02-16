This folder contains a Dockerfile for the Cell Ranger analysis.


## Building the image

In order to build this image, the Cell Ranger software must be downloaded separately to comply with licensing, and should be placed in this folder (`docker/cellranger`).
- The current version of Cell Ranger used here is `9.0.1` and can be downloaded from [10X Genomics Website](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions).
  - Download the `cellranger-9.0.1.tar.gz` file and extract it using the `tar -xf cellranger-9.0.1.tar.gz` command.
  - The extracted directory should be named `cellranger-9.0.1` and contain the Cell Ranger software (we decompress outside the Docker image to save space)

Following download of Cell anger, you can build the image running the following command from this `docker/cellranger` working directory:

```
docker buildx build . -t scpca-cellranger:9.0.1 --platform linux/amd64
```

At this point, the image should be ready for use on the local machine.

## Using the AWS ECR

To use the AWS ECR docker repository, you will need to login  with:

```
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com
```

Note that for the code above to work, you must have set up `aws` command line tools on your machine, and have run [`aws configure`](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html) to set up your credentials for AWS access.
This is a private repository; access is only available to Alex's Lemonade (Data Lab) users at this time.

### Pulling the current image

To pull the current version of the container (skipping the build steps), you can use the following command (once you are logged in):

```
docker pull ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```

### Pushing to AWS ECR

If the updated image needs to be pushed to AWS ECR, you can follow the outline steps below (updating version numbers as needed).
*Do not push an image unless you are sure it is working!*

The current image was pushed with the following commands.

```
docker tag scpca-cellranger:9.0.1 ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
docker tag scpca-cellranger:9.0.1 ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:9.0.1
docker push ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:9.0.1
docker push ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/scpca-cellranger:latest
```
