from django.db import models

class Search(models.Model):
	query = models.CharField(max_length = 200);

	def __str__(self):
		return self.query

class Result(models.Model):
	index = models.IntegerField();
	search = models.ForeignKey(Search, on_delete = models.CASCADE)

	def __str__(self):
		return self.index;