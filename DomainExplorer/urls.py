"""DomainExplorer URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

from DomainExplorer import views

urlpatterns = [
    # Static pages for all the apps
    path('', views.HomePageView.as_view(), name="home"),
    path('download/', views.DownloadPageView.as_view(), name="download-page"),
    path('documentation/', views.DocPageView.as_view(), name="doc-page"),
    path('about/', views.AboutPageView.as_view(), name="about-page"),

    # Include the DIGGER app
    path('', include("domain.urls")),

    path('dummy', views.TemplateView.as_view(template_name='base.html')),
]

# Serve media files during debugging through the django server as well
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)




