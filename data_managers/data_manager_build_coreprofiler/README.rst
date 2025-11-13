This tool downloads and builds the **CoreProfiler** scheme.

You can find the list of available scheme, as well as the reference platforms supported by CoreProfiler, in the <a href="https://gitlab.com/ifb-elixirfr/abromics/coreprofiler/-/blob/main/README.md?ref_type=heads#basic-usage" target="_blank">CoreProfiler documentation</a>.

Please refer to this page for details on how to use the tool and which schema options are available.

Use Galaxy's data manager framework to download and install new CoreProfiler schemes :

If you want to use a scheme from enterobase, you do not need to provide any token or secret.

However, if you want to use a scheme from **pubMLST** or **BigsDB**, you will need to follow a procedure before launching the data manager.

BIGSdb and PubMLST platforms require OAuth1 authentication to access and download the most up-to-date schemes.
While authentication is not strictly mandatory, skipping it may result in downloading outdated schemes.

This authentication involves two types of tokens:

 * Consumer tokens : permanent tokens used to initiate the authentication flow.

 * Access tokens : tokens required to download a scheme.

Procedure for pubMLST schemes (example with borrelia_3-cgMLST-639-pubmlst) :

1. Create an account on the <a href="https://pubmlst.org/bigsdb" target="_blank">pubMLST</a> website.
2. Generate a consumer token and secret from your account settings (My account > API keys > Enter key name > Submit).
3. On your account page, go to Database registrations, check all databases and register.
4. Download <a href="https://gitlab.com/ifb-elixirfr/abromics/coreprofiler" target="_blank">coreprofiler</a> locally and run the following command to obtain your access token and secret

Replace the placeholders with your scheme of interest (example : borrelia_3) your actual consumer token and secret :

"""
coreprofiler db get_request_tokens --scheme <SCHEME_NAME> --consumer_key <YOUR_CONSUMER_TOKEN> --consumer_secret <YOUR_CONSUMER_SECRET>
"""

This command will provide you with a URL to visit in order to authorize client software to access your account. 

After authorizing, it will give you a verification code that you need to enter in the command line prompt.

It will then return your access token and secret.

5. Provide the consumer token, consumer secret, access token, and access secret in the data manager tool by setting this bash variables in a txt file : 

"""
export COREPROFILER_CONSUMER_TOKEN="<YOUR_CONSUMER_TOKEN>"
export COREPROFILER_CONSUMER_SECRET="<YOUR_CONSUMER_SECRET>"
export COREPROFILER_ACCESS_TOKEN="<YOUR_ACCESS_TOKEN>"
export COREPROFILER_ACCESS_SECRET="<YOUR_ACCESS_SECRET>"
"""

6. Set the path to this txt file in your environment by making an environment variable (example : export COREPROFILER_SECRETS_PATH="/path/to/your/secret_file.txt").

Procedure for BigsDB (example with bordetella_1-cgMLST_genus-1415-BIGSdb):

1. Create an account on the <a href="https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?page=registration" target="_blank">BigsDB</a> website.
2. Ask for a consumer token and secret by sending a mail to bigsdb@pasteur.fr (object "API client key").
3. On your account page, go to Database registrations, check all databases and register.
4. Download <a href="https://gitlab.com/ifb-elixirfr/abromics/coreprofiler" target="_blank">coreprofiler</a> locally and run the following command to obtain your access token and secret

Replace the placeholders with your scheme of interest (example : bordetella_1) your actual consumer token and secret :

"""
coreprofiler db get_request_tokens --scheme <SCHEME_NAME> --consumer_key <YOUR_CONSUMER_TOKEN> --consumer_secret <YOUR_CONSUMER_SECRET>
"""

This command will provide you with a URL to visit in order to authorize client software to access your account. 

After authorizing, it will give you a verification code that you need to enter in the command line prompt.

It will then return your access token and secret.

5. Provide the consumer token, consumer secret, access token, and access secret in the data manager tool by setting this bash variables in a txt file : 

"""
export COREPROFILER_CONSUMER_TOKEN="<YOUR_CONSUMER_TOKEN>"
export COREPROFILER_CONSUMER_SECRET="<YOUR_CONSUMER_SECRET>"
export COREPROFILER_ACCESS_TOKEN="<YOUR_ACCESS_TOKEN>"
export COREPROFILER_ACCESS_SECRET="<YOUR_ACCESS_SECRET>"
"""

6. Set the path to this txt file in your environment by making an environment variable (example : export COREPROFILER_SECRETS_PATH="/path/to/your/secret_file.txt").