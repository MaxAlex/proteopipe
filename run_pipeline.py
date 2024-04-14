import os, sys
import yaml
from subprocess import call
import boto3
from botocore.exceptions import ClientError, NoCredentialsError

nextflow_definition = os.path.join(os.path.dirname(__file__), 'workflow', 'workflow_definition.nf')
nextflow_config = os.path.join(os.path.dirname(__file__), 'instance_resources', 'workflow.config')


s3 = boto3.resource('s3')
s3client = boto3.client('s3')



def check_aws_cli_is_logged_in():
    try:
        sts_client = boto3.client('sts')
        identity = sts_client.get_caller_identity()
        account_id = identity['Account']
        user_id = identity['UserId']

        print(f"Found CLI session for AWS account: {account_id}")
        print(f"User ID: {user_id}")
    except NoCredentialsError:
        print("No AWS credentials found.")
    except ClientError as e:
        print(f"Error: {e}")


def check_if_s3_file_exists(bucket, s3path):
    try:
        s3.Object(bucket, s3path).load() # type: ignore
    except:
        return False
    return True


def upload_to_s3(bucket, filename, s3path):
    s3.Bucket(bucket).upload_file(filename, s3path) # type: ignore
    assert(check_if_s3_file_exists(bucket, s3path))


# Note that this will not work on files above a certain size, since
# eventually S3 stops calculating MD5 in a replicable fashion (dependent
# on upload chunks.)
from hashlib import md5
def upload_if_changed(bucket, filename, s3path):
    try:
        headobj = s3client.head_object(Bucket=bucket, Key=s3path)
        s3md5 = headobj['ETag'].strip('"')
    except ClientError: 
        s3md5 = "File not found"
    localmd5 = md5(open(filename, 'rb').read()).hexdigest()
    if s3md5 != localmd5:
        print("Updating %s (%s vs %s)" % (filename, s3md5, localmd5))
        upload_to_s3(bucket, filename, s3path)
    else:
        print("No change to %s (%s)" % (filename, localmd5))


if __name__ == '__main__':
    metaparams = sys.argv[1]
    run_name = os.path.basename(metaparams).split('.')[0]

    assert os.path.exists(metaparams), "Metaparameters file does not exist"
    with open(metaparams, 'r') as f:
        metaparams = yaml.load(f, Loader=yaml.FullLoader)

    # with open('secrets.txt', 'r') as inp:
    #     key, secretkey = inp.read().strip().split()
    #     os.environ['AWS_ACCESS_KEY_ID'] = key
    #     os.environ['AWS_SECRET_ACCESS_KEY'] = secretkey
    check_aws_cli_is_logged_in()

    bucket = metaparams['s3_bucket'].replace('s3://', '')

    s3_paths = {}
    for file_key in ['comet_parameters', 'comet_varmods', 'sage_parameters', 'sage_widewindow_parameters']:
        file_path = metaparams[file_key]
        if file_path.startswith('s3://'):
            assert file_path.startswith(metaparams['s3_bucket']), "All files must be in the same specified s3_bucket"
            # TODO have it auto-copy from the other bucket, if it's not in the right one?
            s3_location = file_path.replace(metaparams['s3_bucket'], '')
        else:
            s3_location = os.path.join(metaparams['s3_parameter_dir'], os.path.basename(file_path))
            upload_if_changed(bucket, file_path, s3_location)

        s3_paths[file_key] = os.path.join('s3://' + bucket, s3_location)

    file_args = sum([['--%s' % k, v] for k, v in s3_paths.items()], [])
    print(file_args)

    if len(sys.argv) > 2 and 'resume' in sys.argv[2]:
        res_arg = ['-resume']
        print("Resuming (%s)" % sys.argv)
    else:
        res_arg = []
        print("From scratch (%s)" % sys.argv)

    cmd = (['nextflow', '-C', nextflow_config, 'run', nextflow_definition, #'-with-dag', 'flowchart.html',
            '--run_name', run_name, '--run_parse_regexp', metaparams['run_parse_regexp'], '--isobaric_type', str(metaparams['isobaric_type']),
            '-bucket-dir', os.path.join(metaparams['s3_bucket'], metaparams['s3_work_dir']), '--fasta_file', metaparams['fasta_file'], 
            '--data_glob', os.path.join(metaparams['s3_bucket'], metaparams['s3_raw_glob'])] + file_args + res_arg)
    print(cmd)
    call(cmd) 
    

