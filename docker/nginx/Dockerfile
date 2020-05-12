FROM nginx:1.17.5

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf.template /etc/nginx/conf.d

COPY entrypoint-nginx.sh /
ENTRYPOINT ["/entrypoint-nginx.sh"]
#ENTRYPOINT ["envsubst", "'${NGINX_PUBLISHED_PATH}'", "<", "/etc/nginx/conf.d/default.conf.template", ">", "/etc/nginx/conf.d/default.conf"]
CMD ["nginx", "-g", "daemon off;"]