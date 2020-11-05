## Set up user credentials on Galaxy to connect to other omero instance

To enable users to set their credentials for this tool,
make sure the file `config/user_preferences_extra.yml` has the following section:

```
    omero_account:
        description: Your OMERO instance connection credentials
        inputs:
            - name: username
              label: Username
              type: text
              required: False
            - name: password
              label: Password
              type:  password
              required: False
```
