# coding=UTF-8
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time,os
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from pyvirtualdisplay import Display

import sys
reload(sys)

sys.setdefaultencoding('utf-8')

display = Display(visible=0, size=(800, 600))
display.start()

username = "clouds"
password = "passwd4BYR"
#title = "post_title is here"
#title = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 
#title = u"测试"
#content = "这是第一行\n\n 这是第二行 \n  这是第三行"
content = "up"

#driver = webdriver.PhantomJS()
driver = webdriver.Firefox()
driver.implicitly_wait(60)
driver.set_page_load_timeout(60)

try:
  driver.get('https://bbs.byr.cn/#!article/Advertising/1905632')
#  driver.get("https://bbs.byr.cn/#!article/test/post")
#  driver.get("https://www.baidu.com")
  driver.execute_script("window.scrollTo(0,document.body.scrollHeight)")
  print " Title : "+driver.title
  #print driver.page_source

  inputElement = driver.find_element_by_name("id")
  inputElement.send_keys(username)
  time.sleep(1.003)
  inputElement = driver.find_element_by_name("passwd")
  inputElement.send_keys(password)
  time.sleep(0.584)
  inputElement.send_keys(Keys.RETURN)

#  print " Login as : "+driver.find_element_by_css_selector("p.username.small").text
  print " Login as : "+driver.find_element_by_css_selector("div.u-login-id a ").text

  time.sleep(3.616)
#  title_input = driver.find_element_by_css_selector("div.input-wrapper.title-input input")
# inputTitle = driver.find_element_by_name("subject")
#  inputTitle.send_keys(title)
#  time.sleep(2.18)
  inputContent = driver.find_element_by_name("content")
  inputContent.click()
  inputContent.send_keys(content)
  time.sleep(3.642)
  print inputContent.get_attribute("value")
  inputContent.send_keys(Keys.CONTROL,Keys.RETURN)
  time.sleep(20.567)
#  ptitle = driver.find_element_by_css_selector("header h3")
#  print " Post title : "+title
  
  #time.sleep(8)

except Exception as e:
  print e
  
finally:
  driver.close()
  display.stop()
  os.system("pkill geckodriver")
  os.system("pkill Xvfb")
