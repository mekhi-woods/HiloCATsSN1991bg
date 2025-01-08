# import os
# import io
# import re
# import sys
# import time
# import requests
# import numpy as np
# import pandas as pd
# import urllib.request
# from requests.auth import HTTPBasicAuth
#
# # BASEURL = "https://fallingstar-data.com/forcedphot"
# # token = "9f57f61e510459775763afe88a2ebfdf2aa37163"
# # headers = {'Authorization': f'Token {token}', 'Accept': 'application/json'}
# # task_url = "https://fallingstar-data.com/forcedphot/queue/1851119/"
#
# API_TOKEN = "9f57f61e510459775763afe88a2ebfdf2aa37163"
#
#
# # pickle = np.load('tmp.npz', allow_pickle=True)
# # data = pickle['data']
#
# data = requests.post('https://star.pst.qub.ac.uk/sne/atlas4/api/objectlist/',
#                      headers={'Authorization': f'Token {API_TOKEN}'},
#                      data={'objectlistid':2}).json()
#
# print(data)
#
# # np.savez('tmp.npz', data=data)
# #
# # count = 0
# # for d in data:
# #     if d['observation_status'] is not None and d['observation_status'].startswith('SN Ia'):
# #         print(d['atlas_designation'],d['observation_status'].replace(' ',''),d['ra'],d['dec'])
# #         count += 1
# #
# #         ids = d['id']
# #         base_url = 'https://star.pst.qub.ac.uk/sne/atlas4/lightcurveforced/1161048951013729300/'
# #         new_url = base_url.replace('1161048951013729300/',str(ids))
# #         print(new_url)
# #
# #         idfile = '/content/ATLAS_ids/' + str(ids)+'.txt'
# #         if os.path.exists(idfile):
# #             continue
# #         urllib.request.urlretrieve(str(new_url), str(idfile))
# #         print(idfile)
# #
# #     if count > 3000:
# #         break
#
#
# # """ Query """
# # task_url = None
# # while not task_url:
# #     with requests.Session() as s:
# #         resp = s.post(f"{BASEURL}/queue/", headers=headers, data={'ra': 44, 'dec': 22, 'mjd_min': 59248.})
# #
# #         if resp.status_code == 201:  # successfully queued
# #             task_url = resp.json()['url']
# #             print(f'The task URL is {task_url}')
# #         elif resp.status_code == 429:  # throttled
# #             message = resp.json()["detail"]
# #             print(f'{resp.status_code} {message}')
# #             t_sec = re.findall(r'available in (\d+) seconds', message)
# #             t_min = re.findall(r'available in (\d+) minutes', message)
# #             if t_sec:
# #                 waittime = int(t_sec[0])
# #             elif t_min:
# #                 waittime = int(t_min[0]) * 60
# #             else:
# #                 waittime = 10
# #             print(f'Waiting {waittime} seconds')
# #             time.sleep(waittime)
# #         else:
# #             print(f'ERROR {resp.status_code}')
# #             print(resp.json())
# #             sys.exit()
# #
# # """ Check Query """
# # result_url = None
# # while not result_url:
# #     with requests.Session() as s:
# #         resp = s.get(task_url, headers=headers)
# #
# #         if resp.status_code == 200:  # HTTP OK
# #             if resp.json()['finishtimestamp']:
# #                 result_url = resp.json()['result_url']
# #                 print(f"Task is complete with results available at {result_url}")
# #                 break
# #             elif resp.json()['starttimestamp']:
# #                 print(f"Task is running (started at {resp.json()['starttimestamp']})")
# #             else:
# #                 print("Waiting for job to start. Checking again in 10 seconds...")
# #             time.sleep(10)
# #         else:
# #             print(f'ERROR {resp.status_code}')
# #             print(resp.json())
# #             sys.exit()
# #
# # with requests.Session() as s:
# #     textdata = s.get(result_url, headers=headers).text
# #
# #     # if we'll be making a lot of requests, keep the web queue from being
# #     # cluttered (and reduce server storage usage) by sending a delete operation
# #     # s.delete(task_url, headers=headers).json()
# #
# # # """ Parse """
# # # dfresult = pd.read_csv(io.StringIO(textdata.replace("###", "")), delim_whitespace=True)
# # # print(dfresult)

# import numpy as np
# import requests
# from requests.auth import HTTPBasicAuth
# import urllib.request
# import os

# API_TOKEN = "7f4e1dee8f52cf0c8cdaf55bf29d04bef4810fb4"
# API_TOKEN = "9f57f61e510459775763afe88a2ebfdf2aa37163"
#
# def main():
#
#         pickle = np.load('tmp.npz', allow_pickle=True)
#         data = pickle['data']
#
#         # data = requests.post(
#         #     'https://star.pst.qub.ac.uk/sne/atlas4/api/objectlist/',
#         #     headers={'Authorization': f'Token {API_TOKEN}'},
#         #     data={'objectlistid':2}
#         # ).json()
#
#         # np.savez('tmp.npz', data=data)
#
#         count = 0
#         for d in data:
#             if d['observation_status'] is not None and d['observation_status'].startswith('SN Ia'):
#                 print(d['atlas_designation'],d['observation_status'].replace(' ',''),d['ra'],d['dec'])
#                 count += 1
#
#                 ids = d['id']
#                 base_url = 'https://star.pst.qub.ac.uk/sne/atlas4/lightcurveforced/1161048951013729300/'
#                 new_url = base_url.replace('1161048951013729300/',str(ids))
#                 print(new_url)
#
#                 idfile = 'ATLAS_ids/' + str(ids)+'.txt'
#                 if os.path.exists(idfile):
#                     continue
#                 urllib.request.urlretrieve(str(new_url), str(idfile))
#                 print(idfile)
#
#             if count > 3000:
#                 break

import requests
import json
from datetime import datetime, timedelta

# Define constants
TNS_API_KEY = "YOUR_TNS_API_KEY"
TNS_BASE_URL = "https://www.wis-tns.org/api/get"  # Base URL for TNS API
# TNS_HEADERS = {
#     "User-Agent": "tns_marker{"
#                   "tns_id": "12345", "
#                   "type": "bot", "
#                   "name": "YourBotName"}"
# }

# Function to build the API request payload
def build_request_payload(params):
    return json.dumps({
        "api_key": TNS_API_KEY,
        "data": params
    })

# Function to query TNS for supernovae discovered by ATLAS within the last year
def query_tns():
    # Define the query parameters
    one_year_ago = (datetime.utcnow() - timedelta(days=365)).strftime("%Y-%m-%d")
    params = {
        "object_type": "supernova",
        "discoverer": "ATLAS",
        "date_start": one_year_ago,
        "classification": "Ia",
        "num_page": 50
    }

    # Send the request
    response = requests.post(
        f"{TNS_BASE_URL}/search",
        headers=TNS_HEADERS,
        data=build_request_payload(params)
    )

    # Handle the response
    if response.status_code == 200:
        data = response.json()
        if data.get("data"):
            return data["data"].get("reply", [])
        else:
            print("No results found.")
            return []
    else:
        print(f"Error: {response.status_code} - {response.text}")
        return []

# Retrieve RA and DEC for each supernova
def get_ra_dec(supernovae):
    ra_dec_list = []
    for sn in supernovae:
        name = sn.get("objname")
        ra = sn.get("radeg")
        dec = sn.get("decdeg")
        if name and ra and dec:
            ra_dec_list.append({"name": name, "ra": ra, "dec": dec})
    return ra_dec_list

# Main function
if __name__ == "__main__":
    supernovae = query_tns()
    if supernovae:
        ra_dec_list = get_ra_dec(supernovae)
        if ra_dec_list:
            print("RA and DEC of Type Ia Supernovae found by ATLAS within the last year:")
            for entry in ra_dec_list:
                print(f"Name: {entry['name']}, RA: {entry['ra']}, DEC: {entry['dec']}")
        else:
            print("No RA and DEC data available.")
    else:
        print("No supernovae found.")