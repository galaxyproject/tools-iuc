# JBrowse2 in Galaxy

    JBrowse is a fast, embeddable genome browser built completely with
    JavaScript and HTML5

Thus, it makes an ideal fit with Galaxy, especially for use as a
workflow summary. E.g. annotate a genome, then visualise all of the
associated datasets as an interactive HTML page. This tool MUST be added
to the "HTML Rendered" allow list in the Galaxy admin panel to function correctly.

## Installation

It is recommended to install this wrapper via the Galaxy Tool Shed.

## Running Locally

The Galaxy tool interface writes out a xml file which is then used to generate
the visualizations. An example used during development/testing can be seen in
`test-data/*/test.xml`. The format is in no way rigorously defined and is
likely to change at any time! Beware. ;)

## Testing with Planemo

The tool comes with standard tests that can be run with `planemo run`.

If you want to make manual tests with `planemo serve`, it is a bit more complex.
If you do it, you will notice that you cannot see any tool output. That's because JBrowse2 needs to be served
by a HTTP server supporting byte-range requests. Planemo does not by default.

Here comes a hacky solution, that works (but is a bit ugly, sorry).

First you need to start a Nginx server that will act as a proxy to the local Galaxy server launched by `planemo serve`. This can be done easily with a docker-compose.yml file like this:

```yaml
services:

  proxy:
    image: quay.io/abretaud/nginx-ldap:latest  # This a basic Nginx image, any other classic Nginx image should work
    volumes:
        - "./nginx/conf:/etc/nginx/conf.d"
        - "/tmp/:/tmp/:ro"
    network_mode: "host"
```

And the mounted `nginx.conf` file:

```text
server {
    listen 80;
    server_name  ~.;

    location /_x_accel_redirect/ {
        internal;
        alias /;
    }

    location / {
        client_max_body_size 50g;
        add_header Accept-Ranges bytes;
        proxy_redirect   http://localhost:9090/       http://$host;
        proxy_pass http://localhost:9090/;
        proxy_force_ranges on;
        proxy_pass_request_headers      on;
    }
}
```

Launch it with `docker compose up -d`.

Then you need to make a small modification of Planemo's code (do it in a dedicated virtualenv), supposing you have installed it in `~/venvs/planemo`, edit `~/venvs/planemo/lib/python*/site-packages/planemo/galaxy/config.py`:

In the `local_galaxy_config` function around line 422, in the `properties.update(...)` call, add this line:

`nginx_x_accel_redirect_base="/_x_accel_redirect",`

Now run `planemo serve`, and then run a `chmod a+rx /tmp/tmpxxxx`, where `/tmp/tmpxxxx` is the directory where Planemo is running Galaxy from.

Now you can access your Planemo Galaxy instance at 2 adresses:

- http://localhost:9090/ : the usual Planemo address, file uploads will work, but you can't access JBrowse2 tool output
- http://localhost:80/ : the same Galaxy server where you *can* access JBrowse2 tool output (don't forget to add the tool in the "HTML Rendered" allow list in the admin panel)

## TODO list

- Support more display options (bug: https://github.com/GMOD/jbrowse-components/issues/5148)
- Reimplement auto coloring of features (custom plugin like https://jbrowse.org/jb2/docs/config_guides/customizing_feature_colors/ ?)
- Support more plugins from https://jbrowse.org/jb2/plugin_store/
- Support Circular view
- Support SV inspector/spreadsheet views
- Investigate using CSI indexing instead of TBI for large chromosomes (https://github.com/GMOD/jbrowse-components/issues/4565#issuecomment-2361153737)
- Synteny: add support for .out (MashMap) .chain (UCSC), .delta (mummer) and .anchors (mcscan) inputs
