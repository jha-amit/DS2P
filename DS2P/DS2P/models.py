from django.db import models

# Create your models here.
class Cost_matrix(models.Model):
    From_node = models.FloatField(blank=True)
    To_node = models.FloatField(blank=True)
    Edge_cost = models.FloatField(blank=True)