import scrapy
from ..items import TestCrawlerItem

# AOS 2020
# Volume 48, Number 5 October 2020:  'https://projecteuclid.org/current/euclid.aos',
# Volume 48, Number 4 (2020), pp. 1875-2503 : 'https://projecteuclid.org/euclid.aos/1597370646'
# Volume 48, Number 3 (2020), pp. 1231-1874: 'https://projecteuclid.org/euclid.aos/1594972813'
# Volume 48, Number 2 (2020), pp. 629-1229: 'https://projecteuclid.org/euclid.aos/1590480025'
# Volume 48, Number 1 Febreuary 2020, pp. 1-627: 'https://projecteuclid.org/euclid.aos/1581930120'


class CoauthorshipSpider(scrapy.Spider):
    name = 'test'
    # start_urls = [  # 2020
    #    'https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-6',  # dec
    #    'https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-5',  # oct
    #    'https://projecteuclid.org/euclid.aos/1597370646',  # aug
    #    'https://projecteuclid.org/euclid.aos/1594972813',  # jun
    #    'https://projecteuclid.org/euclid.aos/1590480025',  # apr
    #    'https://projecteuclid.org/euclid.aos/1581930120'  # feb
    # ]
    start_urls = [  # 2021
        'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-1',  # feb
        # 'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-2', # apr
        # 'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-3', #jun
        # 'https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-4'] # aug
    ]

    def parse(self, response):
        items = TestCrawlerItem()
        all_papers = response.css('.TOCLineItemTextWidth')
        for paper in all_papers:
            summary = paper.css(".TOCLineItemText3::text")[0].extract()
            # print(summary)
            # print(len(summary))
            if len(summary) > 42:  # to skip lines do not contain date
                author = paper.css(".linkBlueToBlue::text").extract()
                # print(author)
                summary_splited = summary.split(',')
                mt_yr = summary_splited[-1]
                mt_yr = mt_yr[mt_yr.find('(')+1:mt_yr.find(')')]
                mt_yr_splited = mt_yr.split(' ')
                month = mt_yr_splited[0]
                year = mt_yr_splited[1]
                print(author, month, year)
                items['author'] = author
                items['month'] = month
                items['year'] = year
        yield items


# extract titles only
# response.css('.TOCLineItemText1::text').extract()
#

# extract names only
# response.css('.TOCLineItemText2::text, .linkBlueToBlue::text').extract)(
#        )
#
# or#
# all_papers=#response.css('.TOCLineItemText2 .linkBlueToBlue')
# for art in all_papers:
#    print(art.css("::text").extract())
#

#
# extract paper, name and time paper by paper#
# all_papers=response.css('.TOCLineItemTextWidth')
# for art in all_papers:
#   author=art.css(".linkBlueToBlue::text").extract()
#   info=art.css(".TOCLineItemText3::text").extract()
#   date=info[0]
#   # keyword = info[1]
#   #print(date)#
#   print(author)
