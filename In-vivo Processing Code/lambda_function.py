
import boto3, csv, json, uuid, os, sys
from datetime import datetime
import numpy as np
from requests_aws4auth import AWS4Auth

from ProcessBMDL import *
from sys import exit
from urllib import request


BUCKET_NAME = os.environ['BUCKET_NAME']
s3_client = boto3.client('s3')

dynamodb = boto3.resource('dynamodb')
counter_table  = dynamodb.Table('Counter_Table')

response = counter_table.put_item(
    Item={  'CounterId': 'todo_counter',
            'counter_value': 0})

def lastFileinDirectory(response):
    file_list = []
    for content in response.get('Contents', []):
        if '.csv' in content['Key'].lower():
            file_list.append(content['Key'])
    filenums = []
    for file in file_list:
        filenum = file.split('.')[0]
        filenums.append(int(filenum))
    filenums.sort()
    return filenums[-1]



def lambda_handler(event, context):

    class GraphqlClient:

        def __init__(self, endpoint, headers):
            self.endpoint = endpoint
            self.headers = headers

        @staticmethod
        def serialization_helper(o):
            if isinstance(o, datetime):
                return o.strftime('%Y-%m-%dT%H:%M:%S.000Z')

        def execute(self, query, operation_name, variables={}):
            data = json.dumps({
                "query": query,
                "variables": variables,
                "operationName": operation_name
            },
                default=self.serialization_helper,
            )
            r = request.Request(
                headers=self.headers,
                url=self.endpoint,
                method='POST',
                data=data.encode('utf8')
            )
            response = request.urlopen(r).read()
            return response.decode('utf8')

    gq_client = GraphqlClient(
        endpoint='https://sav4o2b7vra63eqfxafzktphpu.appsync-api.us-east-2.amazonaws.com/graphql',
        headers={'x-api-key': ' da2-t75j7zwzfbewxjaywrkz4avz6i'}
    )

    response = s3_client.list_objects_v2(Bucket=BUCKET_NAME)
    LAST_FILE_NUM = lastFileinDirectory(response)

    # For each unprocessed csv file, process the file and return output
    csv_filename = str(LAST_FILE_NUM).lower() + '.csv'
    csv_file = s3_client.get_object(Bucket=BUCKET_NAME, Key=csv_filename)
    records_list = csv_file['Body'].read().decode('utf-8-sig').splitlines()   # make a list
    
        
    reverse = False
    if (records_list[0].split(','))[0] == 'CO Concentration': # if reverse order
            reverse = True
    np_records = np.zeros((len(records_list) , 2))
    for idx,element in enumerate(records_list):
        if idx == 0:
            continue
        if reverse:
            np_records[idx,0] = float(element.split(',')[1])//1000 # timestamp
            np_records[idx,1] = float( element.split(',')[0] )*100 /32 # CO Concen
        else:
            np_records[idx,0] = float(element.split(',')[0])//1000 # timestamp
            np_records[idx,1] = float( element.split(',')[1] )*100 /32 # CO Concen

    SUCCESS, volume, r, C, coordinates, refined_timestmaps = ProcessData(csv_filename,np_records)
    if not SUCCESS:
        return None
    volume = str(volume) + ' mL'
    print('...')
    print(f"$$$ VOLUME = {volume}\n")
    print('...')
    print(f"$$$ Rad = {r} mm \n")
    print(f"$$$ Center = {C} mm \n")
    print(f"$$$ Coords = {coordinates} mm \n")
    print(f"$$$ Selected timestamps = {refined_timestmaps} mm \n")
    print('...')
    print('...')
    print('...')


    response = counter_table.update_item(
    Key={'CounterId': 'todo_counter'},
    UpdateExpression='SET counter_value = counter_value + :increment',
    ExpressionAttributeValues={':increment': 1},
    ReturnValues='UPDATED_NEW')
    counter = response['Attributes']['counter_value']
         
    result = gq_client.execute(query= """
            mutation createTodo($createtodoinput: CreateTodoInput!) {
            createTodo(input: $createtodoinput) {
              id
              piezo_no
              bladder_volume
              anterior_d1
              anterior_d2
              anterior_d3
              anterior_d4
              post_d1
              post_d2
              post_d3
              post_d4
              createdAt
              updatedAt
                }
                }
               """,
               
            
        operation_name='createTodo',
        variables={
                    "createtodoinput":
                        {           #id, counter ı dusun
                            "id":str(counter),
                            "piezo_no":0,
                            "bladder_volume":str(volume),
                            "anterior_d1":1,
                            "anterior_d2":2,
                            "anterior_d3":3,
                            "anterior_d4":4,
                            "post_d1":1,
                            "post_d2":2,
                            "post_d3":3,
                            "post_d4":4
                                
                        }
                })


    return None
    

