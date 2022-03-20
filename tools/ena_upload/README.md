# ena-upload-wrapper
Galaxy wrapper for ena-cli-upload

This tool is shipped in a ready to use Galaxy container found [here](https://github.com/ELIXIR-Belgium/ena-upload-container).

## Setting up credentials on Galaxy

The admin of the server can set up global credentials through a file with the format:
```
username: webin_id
password: webin_password
```

The path to this file must be exported as an environment variable called $GALAXY_ENA_SECRETS

Alternatively, the admin can enable users to set their own credentials for this tool.
To enable it, make sure the file `config/user_preferences_extra_conf.yml` has the following section:

```
    ena_webin_account:
        description: Your ENA Webin account details
        inputs:
            - name: webin_id
              label: ENA Webin ID
              type: text
              required: False
            - name: password
              label: Password
              type:  password
              required: False
```
