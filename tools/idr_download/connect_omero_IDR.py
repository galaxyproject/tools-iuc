import requests
from omero.gateway import BlitzGateway

IDR_URL = "https://idr.openmicroscopy.org"


# function to establish an IDR/OMERO connection
def establish_OMERO_connection(uuid_key, usr, psw, host, port=4064):
    # join an existing session if a session key is given, otherwise log in with user/password
    if uuid_key:
        conn = BlitzGateway(host=host, port=port)
        connected = conn.connect(sUuid=uuid_key)
    else:
        conn = BlitzGateway(usr, psw, host=host, port=port)
        connected = conn.connect()
    if not connected:
        raise ConnectionError(f"Could not connect to OMERO server {host}:{port}")
    conn.SERVICE_OPTS.setOmeroGroup('-1')
    return conn


def establish_IDR_connection():
    session = requests.Session()
    response = session.get(f"{IDR_URL}/webclient/?experimenter=-1")
    response.raise_for_status()
    return session
