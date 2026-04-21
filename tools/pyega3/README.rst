PyEGA3
======

Set up user credentials on Galaxy
---------------------------------

To enable users to set their credentials for this tool, make sure the
file ``config/user_preferences_extra.yml`` has the following section:

::

       ega_account:
           description: Your EGA (European Genome Archive) account
           inputs:
               - name: username
                 label: Username
                 type: text
                 required: False
               - name: password
                 label: Password
                 type:  password
                 required: False
