# -*- coding: utf-8 -*-
"""
Created on Sun Sep 13 21:48:58 2020

@author: amit
"""

from django.urls import path
from . import views
urlpatterns = [
    path('', views.home, name='home'),
    path('shortest_path1',views.shortest_path1,name='shortest_path1')
]