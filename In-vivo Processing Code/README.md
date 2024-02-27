# Notes for the processing code:

Runs on the AWS-lambda service, coupled with the S3 service.
Naming convention for the raw data in CSV format are integers in the 
ascending order. Therefore, by selecting the CSV file with the highest
integer name, lambda_function file processes the last CSV uploaded to
the S3 bucket service.

Used for the **Fig. 4e**
