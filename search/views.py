from django.shortcuts import render
from django.http import HttpResponse
# Create your views here.

def query(request) :
	return render(request, 'search/index.html')
	# if(request.method == "POST"):
	# 	return HttpResponse("Received query: " + request.body)
	# return HttpResponse("hello\n");