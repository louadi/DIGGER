from django.contrib import admin
from .models import Gene, Domain, NeaseSavedRun

# Register your models here.
myModels = [Gene, Domain, NeaseSavedRun]  # iterable list
admin.site.register(myModels)