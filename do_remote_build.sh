#!/bin/bash

set -e

INSTANCE_ID="i-0416ac8e7273c149a"

# Start the instance
aws ec2 start-instances --instance-ids $INSTANCE_ID

# Wait for the instance to be running
aws ec2 wait instance-running --instance-ids $INSTANCE_ID

sleep 10 

# Get the public IP address
PUBLIC_IP=$(aws ec2 describe-instances --instance-ids $INSTANCE_ID --query "Reservations[*].Instances[*].PublicIpAddress" --output text)

echo "Public IP Address: $PUBLIC_IP"

rsync -avz -e "ssh -i ~/.ssh/devbox_key.pem -o StrictHostKeyChecking=no" build_resources/ ubuntu@$PUBLIC_IP:~/build_resources

ssh -i ~/.ssh/devbox_key.pem -o StrictHostKeychecking=no ubuntu@$PUBLIC_IP 'cd ~/build_resources && ./docker_image_process.sh | sed "/^Step/d; /^Running/d"' 

if [ "$1" != "keep" ]; then
  # Stop the instance
  aws ec2 stop-instances --instance-ids $INSTANCE_ID
fi
