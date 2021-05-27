#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import requests
import subprocess
from bs4 import BeautifulSoup
#armadillo
req = requests.get("https://a-z-animals.com/animals/armadillo/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-30623")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download
#chicken
req = requests.get("https://a-z-animals.com/animals/chicken/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-47193")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#Mouse
req = requests.get("https://a-z-animals.com/animals/mouse/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-56788")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#Lizard
req = requests.get("https://a-z-animals.com/animals/lizard/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-55458")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#budgerigar
req = requests.get("https://a-z-animals.com/animals/budgerigar/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-46100")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#cat
req = requests.get("https://a-z-animals.com/animals/maine-coon/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-58473")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#cat
req = requests.get("https://a-z-animals.com/animals/anatolian-shepherd-dog/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-20783")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#pig
req = requests.get("https://a-z-animals.com/animals/pig/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img',class_="wp-image-46200")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download

#rabbit
req = requests.get("https://a-z-animals.com/animals/rabbit/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-36734")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download    

#squirrel
req = requests.get("https://a-z-animals.com/animals/squirrel/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-36734")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download        

#bonobo
req = requests.get("https://a-z-animals.com/animals/bonobo/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-45658")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download   

#turkey
req = requests.get("https://a-z-animals.com/animals/turkey/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-46897")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download   

#horse
req = requests.get("https://a-z-animals.com/animals/horse/")
soup = BeautifulSoup(req.content, 'html.parser')
picture = soup.find_all('img', class_="wp-image-39208")
for pictures in picture:
    pictures_url = pictures['src'] # get the href from the tag
    pictures_jpg = [ 'wget', pictures_url ] # just download it using wget.
    subprocess.Popen(pictures_jpg) # run the command to download   




