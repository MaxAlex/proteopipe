import os, sys
import yaml


template_files = [os.path.join('templates', 'workflow.config.template'),
                  os.path.join('templates', 'docker_image_process.sh.template')]

if os.path.exists(sys.argv[1]):
    if sys.argv[1].endswith('.yaml'):
        with open(sys.argv[1], 'r') as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
    else:
        print('Please provide a valid yaml file (or omit the argument to use the interactive mode)')
        sys.exit(1)

    if not all([key in config for key in ['AWS_ACCESS_KEY', 'AWS_SECRET_KEY', 'AWS_REGION', 'ECR_URL', 'ECR_REPO_NAME', 'BATCH_QUEUE_NAME']]):
        print('The provided yaml file does not contain all the necessary keys')
        print('Please provide a yaml file with the following keys: AWS_ACCESS_KEY, AWS_SECRET_KEY, AWS_REGION, ECR_URL, ECR_REPO, BATCH_QUEUE_NAME')
        sys.exit(1)
else:
    config = {}
    config['AWS_ACCESS_KEY'] = input('Enter the AWS access key: ')
    config['AWS_SECRET_KEY'] = input('Enter the AWS secret key: ')
    config['AWS_REGION'] = input('Enter the AWS region: ')
    config['ECR_URL'] = input('Enter the URL for your ECR instance: ')
    config['ECR_REPO_NAME'] = input('Enter the name of the pipeline ECR repository: ')
    config['BATCH_QUEUE_NAME'] = input('Enter the name of the AWS Batch queue: ')


instance_dir = os.path.join(os.path.dirname(__file__), 'instance_resources')
if not os.path.exists(instance_dir):
    os.makedirs(instance_dir)
for template_file in template_files:
    with open(template_file, 'r') as f:
        template = f.read()
    for key in config:
        template = template.replace('{'+key+'}', config[key])
    with open(os.path.join(instance_dir, os.path.basename(template_file).replace('.template', '')), 'w') as f:
        f.write(template)

print('Pipeline configuration files have been created in the instance_resources directory')




