This tool downloads and builds the **CoreProfiler** scheme.
-------------------------------------------------------------------------------

You can find the list of available schemes, as well as the reference platforms supported by CoreProfiler, in the
`CoreProfiler documentation <https://gitlab.com/ifb-elixirfr/abromics/coreprofiler/-/blob/main/README.md?ref_type=heads#basic-usage>`_.

Please refer to this page for details on how to use the tool and which schema options are available.

Use Galaxy's data manager framework to download and install new CoreProfiler schemes.

If you want to use a scheme from **EnteroBase**, you do not need to provide any token or secret.

However, if you want to use a scheme from **pubMLST** or **BigsDB**, you will need to follow a procedure before launching the data manager.

BIGSdb and PubMLST platforms require **OAuth1 authentication** to access and download the most up-to-date schemes.
While authentication is not strictly mandatory, skipping it may result in downloading outdated schemes.

This authentication involves two types of tokens:

* **Consumer tokens**: permanent tokens used to initiate the authentication flow.
* **Access tokens**: tokens required to download a scheme.

Procedure for **pubMLST schemes** (example: ``borrelia_3-cgMLST-639-pubmlst``)
-------------------------------------------------------------------------------

1. Create an account on the `pubMLST website <https://pubmlst.org/bigsdb>`_.
2. Generate a consumer token and secret from your account settings  
   (**My account → API keys → Enter key name → Submit**).
3. On your account page, go to **Database registrations**, check all databases, and register.
4. Download `coreprofiler <https://gitlab.com/ifb-elixirfr/abromics/coreprofiler>`_ locally and run the following command to obtain your access token and secret:

   .. code-block:: bash

      coreprofiler db get_request_tokens --scheme <SCHEME_NAME> \
         --consumer_key <YOUR_CONSUMER_TOKEN> \
         --consumer_secret <YOUR_CONSUMER_SECRET>

   Replace the placeholders with your scheme of interest (example: ``borrelia_3``) and your actual consumer token and secret.

   This command will provide you with a URL to visit in order to authorize the client software to access your account.  
   After authorizing, it will give you a verification code that you need to enter in the command line prompt.  
   It will then return your access token and secret.

5. Provide the consumer token, consumer secret, access token, and access secret in the data manager tool  
   by setting these bash variables in a ``.txt`` file:

   .. code-block:: bash

      export COREPROFILER_CONSUMER_TOKEN="<YOUR_CONSUMER_TOKEN>"
      export COREPROFILER_CONSUMER_SECRET="<YOUR_CONSUMER_SECRET>"
      export COREPROFILER_ACCESS_TOKEN="<YOUR_ACCESS_TOKEN>"
      export COREPROFILER_ACCESS_SECRET="<YOUR_ACCESS_SECRET>"

6. Set the path to this ``.txt`` file in your environment by making an environment variable  
   (example: ``export COREPROFILER_SECRETS_PATH="/path/to/your/secret_file.txt"``).


Procedure for **BigsDB schemes** (example: ``bordetella_1-cgMLST_genus-1415-BIGSdb``)
--------------------------------------------------------------------------------------

1. Create an account on the `BigsDB website <https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?page=registration>`_.
2. Ask for a consumer token and secret by sending an email to ``bigsdb@pasteur.fr``  
   (subject: **API client key**).
3. On your account page, go to **Database registrations**, check all databases, and register.
4. Download `coreprofiler <https://gitlab.com/ifb-elixirfr/abromics/coreprofiler>`_ locally and run the following command to obtain your access token and secret:

   .. code-block:: bash

      coreprofiler db get_request_tokens --scheme <SCHEME_NAME> \
         --consumer_key <YOUR_CONSUMER_TOKEN> \
         --consumer_secret <YOUR_CONSUMER_SECRET>

   Replace the placeholders with your scheme of interest (example: ``bordetella_1``) and your actual consumer token and secret.

   This command will provide you with a URL to visit in order to authorize the client software to access your account.  
   After authorizing, it will give you a verification code that you need to enter in the command line prompt.  
   It will then return your access token and secret.

5. Provide the consumer token, consumer secret, access token, and access secret in the data manager tool  
   by setting these bash variables in a ``.txt`` file:

   .. code-block:: bash

      export COREPROFILER_CONSUMER_TOKEN="<YOUR_CONSUMER_TOKEN>"
      export COREPROFILER_CONSUMER_SECRET="<YOUR_CONSUMER_SECRET>"
      export COREPROFILER_ACCESS_TOKEN="<YOUR_ACCESS_TOKEN>"
      export COREPROFILER_ACCESS_SECRET="<YOUR_ACCESS_SECRET>"

6. Set the path to this ``.txt`` file in your environment by making an environment variable  
   (example: ``export COREPROFILER_SECRETS_PATH="/path/to/your/secret_file.txt"``).
