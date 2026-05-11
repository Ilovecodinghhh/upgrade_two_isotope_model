import boto3
import multiprocessing
import os
from botocore import UNSIGNED
from botocore.client import Config
# First, set up access to S3 without credentials.
bucket_name = "blended-tropomi-gosat-methane"
s3 = None
def initialize():
    global s3
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
initialize()
# Loop through the folders for each months and collect all S3 paths
months = ([f"2018-{str(m).zfill(2)}" for m in range(4,13)] +
          [f"{y}-{str(m).zfill(2)}" for m in range(1,13) for y in range(2019,2024)])

s3_paths = []
for month in months:
    Prefix=(f"data/{month}/")
    for key in s3.list_objects(Bucket=bucket_name, Prefix=Prefix)["Contents"]:
        s3_paths.append(key["Key"])
print(f"Going to download {len(s3_paths)} files.")

# Download the files using multiple cores
storage_dir = "/n/holyscratch01/jacob_lab/nbalasus/test"
os.makedirs(storage_dir, exist_ok=True)

def download_from_s3(s3_path):
    file =  os.path.basename(s3_path)
    local_file_path = os.path.join(storage_dir,file)
    s3.download_file(bucket_name, s3_path, local_file_path)

with multiprocessing.Pool(112, initialize) as pool:
    pool.map(download_from_s3, s3_paths)
    pool.close()
    pool.join()