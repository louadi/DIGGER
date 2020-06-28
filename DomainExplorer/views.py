from django.views.generic import TemplateView


class HomePageView(TemplateView):
    template_name = "home.html"


class DownloadPageView(TemplateView):
    template_name = "download.html"


class DocPageView(TemplateView):
    template_name = "documentation.html"


class AboutPageView(TemplateView):
    template_name = "about.html"
