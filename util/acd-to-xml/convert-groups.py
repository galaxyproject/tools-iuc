import fileinput
import json


def getFromDict(data, ml):
    # http://stackoverflow.com/a/14692747
    return reduce(lambda d, k: d[k], ml, data)


def setInDict(data, ml, v):
    getFromDict(data, ml[:-1])[ml[-1]] = v

dsc = {}
for line in fileinput.input():
    if not line.startswith('#'):
        data = line.split()
        edam_topic = data[0]
        category = data[1]
        description = ' '.join(data[2:])

        path = category.split(':')
        setInDict(dsc, path, {
            '_description': description,
            '_edam': edam_topic,
        })

print json.dumps(dsc)
