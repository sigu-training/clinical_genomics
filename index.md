---
layout: index
title: Home
---



## Welcome to the documentation site of the           
## Training Course on [“Data analysis and interpretation for clinical genomics” using Galaxy](https://sigu-training.github.io/clinical_genomics)

___

#### Here you can find a tutorials collection to support the course.

{::nomarkdown}

{% assign pages_list = site.pages | sort:"url" %}
    {% for node in pages_list %}
      {% if node.title != null %}
        {% if node.layout == "page" %}
          <a class="sidebar-nav-item{% if page.url == node.url %} active{% endif %}" href="{{site.url}}{{ node.url }}">{{ node.title }}
          <p class="note">{{node.summary}}</p></a>
        {% endif %}
      {% endif %}
    {% endfor %}
{:/}


	  
