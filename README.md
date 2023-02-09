# Single-Cell Analysis of RNA Browser
This repo contains the code for the SCARAB browser for exploring the data from the manuscript Cheng, J. B., Sedgewick, A. J., Finnegan, A. I., Harirchian, P., Lee, J., Kwon, S., … Cho, R. J. (2018). Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution. Cell Reports, 25(4), 871–883. http://doi.org/10.1016/j.celrep.2018.09.006

`app.R` contains the main version for use in a full sized browser window,
 and `app_narrow.R` has a smaller version for use in an inset window.

## Setup ##
Build image from Dockerfile
Map container directories:
  - apps in /srv/shiny-server/
  - logs in /var/log/shiny-server


```
docker build -f scarab_docker/Dockerfile -t scarab scarab_docker (needs shiny-server.sh)
docker run -d -p 3838:3838 -v /home/ajsedgewick/scarab_vol/apps/:/srv/shiny-server/ -v /home/ajsedgewick/scarab_vol/logs:/var/log/shiny-server/ -v /home/ajsedgewick/scarab_vol/:/home/shiny/ -v /data/ajsedgewick:/data/ registry:5000/ajsedgewick/scarab:latest shiny-server.sh
```

Set up proxy following this tutorial: https://support.rstudio.com/hc/en-us/articles/213733868-Running-Shiny-Server-with-a-Proxy

We're using nginx so set up scarab-research block in /etc/nginx/sites-enabled/research looks like this

```
# scarab-research:443
server {
    listen  443 ;
        server_name scarab-research.nantomics.com;

    ### SSL cert files ###
    ssl on;
    ssl_certificate      /etc/ssl/research.nantomics.com.crt;
    ssl_certificate_key  /etc/ssl/research.nantomics.com.key;
    ssl_dhparam         /etc/ssl/dhparams.pem;
    ssl_protocols       TLSv1 TLSv1.1 TLSv1.2;
    ssl_ciphers         HIGH:!aNULL:!MD5:!3DES;
    ssl_prefer_server_ciphers on;
    keepalive_timeout   60;
    ssl_session_cache   shared:SSL:10m;
    ssl_session_timeout 10m;

    access_log /var/log/nginx/access.log;
    error_log /var/log/nginx/error.log;

    root /var/www/html/scarab;

    index index.html;
    location / {
      # First attempt to serve request as file, then
      # as directory, then fall back to displaying a 404.
      #try_files $uri $uri/ =404;

      #redirect to main so I can develop on scdev
      rewrite ^/$ /main break;

      proxy_set_header    Host $host;
      proxy_set_header    X-Real-IP $remote_addr;
      proxy_set_header    X-Forwarded-For $proxy_add_x_forwarded_for;
      proxy_set_header    X-Forwarded-Proto $scheme;

      proxy_pass          http://localhost:3838;
      proxy_read_timeout  90;
      proxy_buffering off;

      proxy_redirect      http://localhost:3838/ $scheme://$host/;
    }
}
```