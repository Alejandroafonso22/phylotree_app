import pandas as pd
import requests
import subprocess
from bs4 import BeautifulSoup
import base64



specie = "gallus gallus"
req = requests.get("https://www.google.com/search?q="+specie+"&tbm=isch&hl=es&sa=X&ved=2ahUKEwj0mOH6xtPwAhUJdxoKHTFcAgsQBXoECAEQOw&biw=1275&bih=942")
soup = BeautifulSoup(req.content, "html.parser")
picture = soup.find_all('img', src_="")[1]
picture_with_image = picture['src']