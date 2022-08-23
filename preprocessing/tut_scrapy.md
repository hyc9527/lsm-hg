# A companion to the [tutorial](https://www.youtube.com/watch?v=ve_0h4Y8nuI&list=PLhTjy8cBISEqkN-5Ku_kXG4QW33sxQo0t&index=1) on Scrapy

### 1 - Web Scraping, Spiders and Crawling

How do we extract information (like title, authors, price, ...) about [recent best sellers in Amazon Books section](https://www.amazon.com/Best-Sellers-Books/zgbs/books/ref=zg_bs_nav_0)?


Sol: devide it into four tasks:
- [task 0](#task-0):   Crawling titles, authors and quotes from a simple [webpage](https://quotes.toscrape.com/).
- [task 1](#task-1)
- [task 2](#task-2)
- [task 3](#task-3)


--------------------

#### Task 0

###### Objective:

Crawling basic data  from a simple [webpage](https://quotes.toscrape.com/) with help of  python module `Scrapy`.


Tools/concepts:
- scrapy module
    - create a project
    - write a spider script
- html/css path
    - SelectorGadget, a chrome extension
    - css selector
    - xpath selector



| Target        | valid css path |
|---------------|----------------|
| webpage title |                |
| quotes        |                |
| authors       |                |
| 'next' button |                |


###### Steps:
The task is divided into three steps:
1. find valid paths of target information.
2. ceate a scrapy project, including project folders and spider script.
3. run spider.

Now we start with scratch.
Check module `Scrapy`  is installed.

```sh
   python3 -m scrapy
```

Create a project in a shell

```sh
    scrapy startproject my_quote_crawler
```


Create a spider script in the folder `spiders`

```sh
    touch ./my_quote_crawler/my_quote_crawler/spiders/quote_spider.py
```

The spider 'quote_spider.py' is to find title of the webpage, which is  "Quotes to Scrape".

```python
import scrapy
class QuoteSpider(scrapy.Spider):
    name = 'quote' # spider name
    start_urls = ['https://quotes.toscrape.com/']
    def parse(self, response):
        title = response.css('title').extract() # extract title info from the webpage
        yield {'titletext':title}
```

To use the spider in a shell, go to the [project folder](/Users/yichen/Dropbox/s21/research/weekly/week19/mcem-real-data/tut_scrapy/my_quote_crawler) and then run the spider in cmd to check if the title is found.

```sh
cd ./my_quote_crawler
scrapy crawl quote
```

Interactive env

```
 scrapy shell 'https://quotes.toscrape.com/'
```

Using css selector, there is a chrome extension "SelectorGadget" to find css conditions with wihch data is extracted. see [lect 9](https://www.youtube.com/watch?v=FQv-whbCfKs&list=PLhTjy8cBISEqkN-5Ku_kXG4QW33sxQo0t&index=9)
To extract title from the webpage:

```
In : response.css("title::text")
Out: [<Selector xpath='descendant-or-self::title/text()' data='Quotes to Scrape'>]

In : response.css("title")
Out: [<Selector xpath='descendant-or-self::title' data='<title>Quotes to Scrape</title>'>]

In: response.css("title::text").extract()
Out: ['Quotes to Scrape']

In: response.css("title::text")[0].extract()
Out: 'Quotes to Scrape'

```

To extract quotes,

```
In : response.css("span.text::text").extract()
Out:
['“The world as we have created it is a process of our thinking. It cannot be changed without changing our thinking.”',
'“It is our choices, Harry, that show what we truly are, far more than our abilities.”',
... ...
'“A day without sunshine is like, you know, night.”']
```

To extract authors,

```
In [10]: response.css(".author::text").extract()
Out[10]:
['Albert Einstein', 'J.K. Rowling', 'Albert Einstein', 'Jane Austen', 'Marilyn Monroe', 'Albert Einstein','André Gide','Thomas A. Edison','Eleanor Roosevelt','Steve Martin']
```

Using xpath selector, button information can be crawled. see [lect 10](https://www.youtube.com/watch?v=INm8yR4aYjk&list=PLhTjy8cBISEqkN-5Ku_kXG4QW33sxQo0t&index=10)

```
scrapy shell 'https://quotes.toscrape.com/'
response.xpath("//title").extract()
response.xpath("//title/text()").extract()
response.xpath("//span[@class='text']/text()").extract()
response.xpath("//span[@class='text']/text()")[0].extract()
response.css("li.next a").xpath("@href").extract()
response.css("//a").xpath("@href").extract()
```
Note that, in css,  single quote must match with single and double quote with double.



With what we have learned about css selector so far, go back to the spider script 'quote_spider.py'. Make a change to extract quotes and authors.
```python
import scrapy
class QuoteSpider(scrapy.Spider):
    name = 'quote' # spider name
    start_urls = ['https://quotes.toscrape.com/']
    #def parse(self, response):
    #    title = response.css('title').extract() # extract title info from the webpage
    #    yield {'titletext':title}
    def parse(self, response):
        all_div_quotes = response.css('div.quote')[0] # for test
        #all_div_quotes = response.css('div.quote')
        title = all_div_quotes.css('span.text::text').extract()
        author=all_div_quotes.css('.author::text').extract()
        tag = all_div_quotes.css('.tag::text').extract()
        yield{
            'title': title,
            'author':author,
            'tag':tag
        }
```

Run the spider in cmd to check out the result.
```
scrapy crawl quote
```

[lect 13](https://www.youtube.com/watch?v=kkWhQKtxT2I&list=PLhTjy8cBISEqkN-5Ku_kXG4QW33sxQo0t&index=13)
It is more convinient and practical to store extracted data as temporary container, called item. Set it up in the file `items.py`

``` python
# Define here the models for your scraped items
#
# See documentation in:
# https://docs.scrapy.org/en/latest/topics/items.html

import scrapy
class MyQuoteCrawlerItem(scrapy.Item):
    # define the fields for your item here like:
    title = scrapy.Field()
    author = scrapy.Field()
    tag = scrapy.Field()
    pass
```
add a line in `quote_spider.py`
``` python
from ..items import MyQuoteCrawlerItemy

```

Store data in a json/csv/xml file
```
scrapy crawl quote -o output.json
scrapy crawl quote -o output.csv
scrapy crawl quote -o output.xml
```

