import json

with open("meta.json") as jh:
    data = json.load(jh)

print("ID\tcollection_date\tcountry\tstate\tlocality")
for k in data:
    collection_date = data[k]['collected']
    location = data[k]['location']
    print("%s\t%s\t%s\t%s\t%s" % (k, collection_date, location['country'], location['state'], location['locality']))
