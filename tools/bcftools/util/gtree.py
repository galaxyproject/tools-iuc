import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name=__name__)
import xml.etree.ElementTree as etree

starter = {
    'name': -100,
    'label': -99,
    'type': -80,

    'input': -52,
    'format': -51,
    'truevalue': -30,
    'falsevalue': -29,

    # Adjust for lower case ascii ORDs
    'selected': 30 + 97,
    # All the way to the end
    'help': 31 + 97,
}

def priority(key):
    if key in starter:
        return starter[key]
    else:
        return ord(key[0])

class GElementTree(etree.ElementTree):

    #http://stackoverflow.com/questions/2741480/can-elementtree-be-told-to-preserve-the-order-of-attributes
    def _serialize_xml(self, write, elem, encoding, qnames, namespaces):
        if write is None:
            raise Exception("Blah")
        tag = elem.tag
        text = elem.text
        if tag is etree.Comment:
            write("<!--%s-->" % etree._encode(text, encoding))
        elif tag is etree.ProcessingInstruction:
            write("<?%s?>" % etree._encode(text, encoding))
        else:
            tag = qnames[tag]
            if tag is None:
                if text:
                    write(etree._escape_cdata(text, encoding))
                for e in elem:
                    self._serialize_xml(write, e, encoding, qnames, None)
            else:
                write("<" + tag)
                items = elem.items()
                if items or namespaces:
                    if namespaces:
                        for v, k in sorted(namespaces.items(),
                                        key=lambda x: x[1]):  # sort on prefix
                            if k:
                                k = ":" + k
                            write(" xmlns%s=\"%s\"" % (
                                k.encode(encoding),
                                etree._escape_attrib(v, encoding)
                                ))
                    #for k, v in sorted(items):  # lexical order
                    #for k, v in items: # Monkey patch
                    for k, v in sorted(items, key=lambda x: priority(x[0])):
                        if isinstance(k, etree.QName):
                            k = k.text
                        if isinstance(v, etree.QName):
                            v = qnames[v.text]
                        else:
                            v = etree._escape_attrib(v, encoding)
                        write(" %s=\"%s\"" % (qnames[k], v))
                if text or len(elem):
                    write(">")
                    if text:
                        write(etree._escape_cdata(text, encoding))
                    for e in elem:
                        self._serialize_xml(write, e, encoding, qnames, None)
                    write("</" + tag + ">")
                else:
                    write(" />")
        if elem.tail:
            write(etree._escape_cdata(elem.tail, encoding))

    def gwrite(self, file_or_filename):
        if hasattr(file_or_filename, "write"):
            file = file_or_filename
        else:
            file = open(file_or_filename, "wb")

        write = file.write
        write("<?xml version='1.0' encoding='%s'?>\n" % 'utf-8')

        qnames, namespaces = etree._namespaces(
            self._root, 'utf-8', None
        )
        serialize = self._serialize_xml
        serialize(write, self._root, 'utf-8', qnames, namespaces)

        if file_or_filename is not file:
            file.close()
