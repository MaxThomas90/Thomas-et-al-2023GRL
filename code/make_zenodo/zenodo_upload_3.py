#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 19:24:09 2023

@author: maxthomas_UO
"""

import requests
import json
import sys

#suite = sys.argv[1]
#if suite == 'u-ci501':
#    ensemble_member = 1
#elif suite == 'u-cm483':
#    ensemble_member = 2
#elif suite == 'u-cn043':
#    ensemble_member = 3
#elif suite == 'u-cn077':
#    ensemble_member = 4
#description = 'Replication data for: Future response of Antarctic continental shelf temperatures to ice shelf basal melting and calving (SSP585FW ensemble member %s, suite id %s). These files contain atmosphere, ocean, and sea-ice data produced by HadGEM3-GC3.1-LL, subset to provide regional and vertical slices that are sufficient to reproduce the manuscript figures and analysis.' % (ensemble_member,suite)
#title = 'Replication data for: Future response of Antarctic continental shelf temperatures to ice shelf basal melting and calving (SSP585FW ensemble member %s, suite id %s).' % (ensemble_member, suite)

ACCESS_TOKEN = 'e4vgvVosPNiKiKYovVSTlBT39UMQnvOKjWFpsmba0XwfVyiM2upxGJDG8N3d' # main
#ACCESS_TOKEN = 'jVBdAKHNp65OAaoXv83FbTNe59PZiyi8drRr7Hwidb7beaVsVAMVLjaeC5nE' # sandbox

headers = {"Content-Type": "application/json"}
params = {'access_token': ACCESS_TOKEN}
r = requests.post('https://zenodo.org/api/deposit/depositions',
#r = requests.post('https://sandbox.zenodo.org/api/deposit/depositions',
                   params=params,
                   json={},
                   headers=headers)
r.status_code
# 201
zenodo_dict = r.json()
print(zenodo_dict)

bucket_url = zenodo_dict["links"]["bucket"]
print(bucket_url)
''' New API '''
for suite in ['SSP585FW_ensemble']:
    filename = suite + '.tar.gz'
    path = "/nesi/nobackup/nesi00442/thoma97p/zips_for_zenodo/%s" % filename
    print(filename)

    with open(path, "rb") as fp:
        r = requests.put(
            "%s/%s" % (bucket_url, filename),
            data=fp,
            params=params,
        )
        print(suite)
#        print(r.json())


'''
Old API
Get the deposition id from the previous response
'''
deposition_id = zenodo_dict['id']
data = {'name': 'readme'}
files = {'file': open('readme', 'rb')}
r = requests.post('https://zenodo.org/api/deposit/depositions/%s/files' % deposition_id,
                    params={'access_token': ACCESS_TOKEN}, data=data,
                    files=files)
r.status_code
# 201
r.json()

#data = {
#     'metadata': {
#         'title': title,
#         'upload_type': 'replication data',
#         'description': description,
#         'creators': [{'name': 'Thomas, Max',
#                       'affiliation': 'University of Otago'}]
#     }
# }

#r = requests.put('https://zenodo.org/api/deposit/depositions/%s' % deposition_id,
#                  params={'access_token': ACCESS_TOKEN}, data=json.dumps(data),
#                  headers=headers)
#r.status_code

print('DEPOSITION_ID = ' + str(deposition_id))

