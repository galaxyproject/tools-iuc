Webin CLI wrapper
======================
the submission of genome assemblies using the 
[Webin-cli](https://github.com/enasequence/webin-cli)

This is a Galaxy wrapper for `webin cli <https://github.com/enasequence/webin-cli>`__  tool from `ENA <https://www.ebi.ac.uk/ena/browser/home>`__.


Setting up credentials on Galaxy
--------------------------------

Enable users to set their own credentials
for this tool. To enable it, make sure the file
``config/user_preferences_extra_conf.yml`` has the following section:

.. code-block:: yaml

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
