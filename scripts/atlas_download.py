import os
import io
import re
import sys
import time
import requests
import numpy as np
import pandas as pd
import urllib.request
from requests.auth import HTTPBasicAuth

# BASEURL = "https://fallingstar-data.com/forcedphot"
# token = "9f57f61e510459775763afe88a2ebfdf2aa37163"
# headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}
# task_url = "https://fallingstar-data.com/forcedphot/queue/1851119/"

API_TOKEN = "9f57f61e510459775763afe88a2ebfdf2aa37163"


# pickle = np.load('tmp.npz', allow_pickle=True)
# data = pickle['data']

data = requests.post('https://star.pst.qub.ac.uk/sne/atlas4/api/objectlist/',
                     headers={'Authorization': f'Token {API_TOKEN}'},
                     data={'objectlistid':2}).json()

print(data)

# np.savez('tmp.npz', data=data)
#
# count = 0
# for d in data:
#     if d['observation_status'] is not None and d['observation_status'].startswith('SN Ia'):
#         print(d['atlas_designation'],d['observation_status'].replace(' ',''),d['ra'],d['dec'])
#         count += 1
#
#         ids = d['id']
#         base_url = 'https://star.pst.qub.ac.uk/sne/atlas4/lightcurveforced/1161048951013729300/'
#         new_url = base_url.replace('1161048951013729300/',str(ids))
#         print(new_url)
#
#         idfile = '/content/ATLAS_ids/' + str(ids)+'.txt'
#         if os.path.exists(idfile):
#             continue
#         urllib.request.urlretrieve(str(new_url), str(idfile))
#         print(idfile)
#
#     if count > 3000:
#         break


# """ Query """
# task_url = None
# while not task_url:
#     with requests.Session() as s:
#         resp = s.post(f"{BASEURL}/queue/", headers=headers, data={'ra': 44, 'dec': 22, 'mjd_min': 59248.})
#
#         if resp.status_code == 201:  # successfully queued
#             task_url = resp.json()['url']
#             print(f'The task URL is {task_url}')
#         elif resp.status_code == 429:  # throttled
#             message = resp.json()["detail"]
#             print(f'{resp.status_code} {message}')
#             t_sec = re.findall(r'available in (\d+) seconds', message)
#             t_min = re.findall(r'available in (\d+) minutes', message)
#             if t_sec:
#                 waittime = int(t_sec[0])
#             elif t_min:
#                 waittime = int(t_min[0]) * 60
#             else:
#                 waittime = 10
#             print(f'Waiting {waittime} seconds')
#             time.sleep(waittime)
#         else:
#             print(f'ERROR {resp.status_code}')
#             print(resp.json())
#             sys.exit()
#
# """ Check Query """
# result_url = None
# while not result_url:
#     with requests.Session() as s:
#         resp = s.get(task_url, headers=headers)
#
#         if resp.status_code == 200:  # HTTP OK
#             if resp.json()['finishtimestamp']:
#                 result_url = resp.json()['result_url']
#                 print(f"Task is complete with results available at {result_url}")
#                 break
#             elif resp.json()['starttimestamp']:
#                 print(f"Task is running (started at {resp.json()['starttimestamp']})")
#             else:
#                 print("Waiting for job to start. Checking again in 10 seconds...")
#             time.sleep(10)
#         else:
#             print(f'ERROR {resp.status_code}')
#             print(resp.json())
#             sys.exit()
#
# with requests.Session() as s:
#     textdata = s.get(result_url, headers=headers).text
#
#     # if we'll be making a lot of requests, keep the web queue from being
#     # cluttered (and reduce server storage usage) by sending a delete operation
#     # s.delete(task_url, headers=headers).json()
#
# # """ Parse """
# # dfresult = pd.read_csv(io.StringIO(textdata.replace("###", "")), delim_whitespace=True)
# # print(dfresult)