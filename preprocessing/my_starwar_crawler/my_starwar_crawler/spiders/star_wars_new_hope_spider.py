import scrapy
from ..items import mySpiderItem


class mySpider(scrapy.Spider):
    name = 'new_hope'
    start_urls = [  # 2021
      'https://imsdb.com/scripts/Star-Wars-A-New-Hope.html',  # Episode IV A New Hope
    ]

    def parse(self, response):
        items = mySpiderItem()
        script = response.css('b::text')[5:]
        for frame in script:
            #data = frame.extract()[:-2]

            data = frame.extract()
            print(data)
            items['data'] = data
            #print(items)
            yield items
