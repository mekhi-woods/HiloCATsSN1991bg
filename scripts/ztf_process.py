def ztf_collection(submit=False, limit=1000):
    objs = dict_unpacker(ATLAS_SAVE_TXT)
    print("Number of targets =", len(objs))
    ralist, declist, jdslist, jdelist = [], [], [], []
    i = 0
    for obj in objs:
        ralist.append(float(objs[obj]['ra']))
        declist.append(float(objs[obj]['dec']))
        jdslist.append(float(objs[obj]['MJDs'])-100)
        jdelist.append(float(objs[obj]['MJDe'])+100)
        i += 1
        if i % limit == 0:
            if submit:
                ztf_submit_post(ralist, declist, jdslist, jdelist)
            ralist, declist, jdslist, jdelist = [], [], [], []
    if submit and len(ralist) > 0:
        ztf_submit_post(ralist, declist, jdslist, jdelist)

    return

def ztf_submit_post(ra_list, dec_list, jds, jde):
    print('Submiting request to ZTF...')

    email = 'mekhidw@hawaii.edu' # email you subscribed with.
    userpass = 'wxdk286' # password that was issued to you.

    ra, dec = json.dumps(ra_list), json.dumps(dec_list)
    jdend, jdstart = json.dumps(jde), json.dumps(jds)
    payload = {'ra': ra, 'dec': dec, 'jdstart': jdstart, 'jdend': jdend, 'email': email, 'userpass': userpass}

    # fixed IP address/URL where requests are submitted:
    url = 'https://ztfweb.ipac.caltech.edu/cgi-bin/batchfp.py/submit'
    r = requests.post(url, auth=('ztffps', 'dontgocrazy!'), data=payload)

    print('RA ['+str(type(payload['ra'][0]))+']:', payload['ra'])
    print('DEC ['+str(type(payload['dec'][0]))+']:', payload['dec'])
    print('MJD Start ['+str(type(payload['jdstart'][0]))+']:', payload['jdstart'])
    print('MJD End ['+str(type(payload['jdend'][0]))+']:', payload['jdend'])
    print("Status_code=", r.status_code)
    return

def ztf_alt_collection():
    ra, dec, jds, jde = 186.07860833333334, 10.446222222222222, 2460350, 2460450

    email = 'mekhidw@hawaii.edu' # email you subscribed with.
    userpass = 'wxdk286' # password that was issued to you.

    cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={dec}&dec={ra}&jdstart={jds}&jdend={jde}&email={email}&userpass={userpass}\""
    print(cmd)
    os.system(cmd)

    return
