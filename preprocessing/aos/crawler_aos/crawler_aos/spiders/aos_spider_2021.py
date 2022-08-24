import scrapy
from ..items import TestCrawlerItem


class CoauthorshipSpider(scrapy.Spider):
    name = 'aos-2021'
    #start_urls = [  # 2020
    #    'https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-6',  # dec
    #    'https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-5',  # oct
    #    'https://projecteuclid.org/euclid.aos/1597370646',  # aug
    #    'https://projecteuclid.org/euclid.aos/1594972813',  # jun
    #    'https://projecteuclid.org/euclid.aos/1590480025',  # apr
    #    'https://projecteuclid.org/euclid.aos/1581930120'  # feb
    # ]

    start_urls = [  # 2021
        'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-1',  # feb
        'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-2', # apr
        'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-3', #jun
        'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-4' # aug
    ]

    def parse(self, response):
        items = TestCrawlerItem()
        all_papers = response.css('.TOCLineItemTextWidth')
        for paper in all_papers:
            summary = paper.css(".TOCLineItemText3::text")[0].extract()
            if len(summary) > 42:  # to skip lines do not contain date and authors' names
                summary_splited = summary.split(',')
                mt_yr = summary_splited[-1]
                mt_yr = mt_yr[mt_yr.find('(')+1:mt_yr.find(')')]
                mt_yr_splited = mt_yr.split(' ')
                month = mt_yr_splited[0]
                year = mt_yr_splited[1]
                author = paper.css(".linkBlueToBlue::text").extract()
                items['author'] = author.copy()
                items['month'] = month
                items['year'] = year
                # print(items)
                yield items






