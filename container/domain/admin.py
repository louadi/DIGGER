from django.contrib import admin
from .models import Gene, Domain, NeaseSavedRun

class NeaseSavedRunAdmin(admin.ModelAdmin):
    readonly_fields = ["run_id"]
    list_display = ["run_id", "timestamp_of_creation", "saved_for_days", "organism", "input_format"]
    search_fields = ['run_id']
    list_filter = ["organism", "input_format", "saved_for_days"]
    ordering = ["-timestamp_of_creation"]

# Register your models here.
myModels = [Gene, Domain]  # iterable list
admin.site.register(myModels)
admin.site.register(NeaseSavedRun, NeaseSavedRunAdmin)