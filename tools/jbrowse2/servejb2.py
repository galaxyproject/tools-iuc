#!/usr/bin/env python3

# spec: simplest python web server with range support and multithreading that takes root path,
# port and bind address as command line arguments; by default uses the current dir as webroot,
# port 8000 and bind address of 0.0.0.0
# borrowed from https://github.com/danvk/RangeHTTPServer
# and reborrowed from https://gist.github.com/glowinthedark/b99900abe935e4ab4857314d647a9068


import argparse
import functools
import os
import re
import socketserver
from http.server import SimpleHTTPRequestHandler
import webbrowser

DEFAULT_PORT = 8080


def copy_byte_range(infile, outfile, start=None, stop=None, bufsize=16 * 1024):
    """Like shutil.copyfileobj, but only copy a range of the streams.

    Both start and stop are inclusive.
    """
    if start is not None:
        infile.seek(start)
    while 1:
        to_read = min(bufsize, stop + 1 - infile.tell() if stop else bufsize)
        buf = infile.read(to_read)
        if not buf:
            break
        outfile.write(buf)


BYTE_RANGE_RE = re.compile(r"bytes=(\d+)-(\d+)?$")


def parse_byte_range(byte_range):
    """Returns the two numbers in 'bytes=123-456' or throws ValueError.

    The last number or both numbers may be None.
    """
    if byte_range.strip() == "":
        return None, None

    m = BYTE_RANGE_RE.match(byte_range)
    if not m:
        raise ValueError("Invalid byte range %s" % byte_range)

    first, last = [x and int(x) for x in m.groups()]
    if last and last < first:
        raise ValueError("Invalid byte range %s" % byte_range)
    return first, last


class RangeRequestHandler(SimpleHTTPRequestHandler):
    """Adds support for HTTP 'Range' requests to SimpleHTTPRequestHandler

    The approach is to:
    - Override send_head to look for 'Range' and respond appropriately.
    - Override copyfile to only transmit a range when requested.
    """

    def handle(self):
        try:
            SimpleHTTPRequestHandler.handle(self)
        except Exception:
            # ignored, thrown whenever the client aborts streaming (broken pipe)
            pass

    def send_head(self):
        if "Range" not in self.headers:
            self.range = None
            return SimpleHTTPRequestHandler.send_head(self)
        try:
            self.range = parse_byte_range(self.headers["Range"])
        except ValueError:
            self.send_error(400, "Invalid byte range")
            return None
        first, last = self.range

        # Mirroring SimpleHTTPServer.py here
        path = self.translate_path(self.path)
        f = None
        ctype = self.guess_type(path)
        try:
            f = open(path, "rb")
        except IOError:
            self.send_error(404, "File not found")
            return None

        fs = os.fstat(f.fileno())
        file_len = fs[6]
        if first >= file_len:
            self.send_error(416, "Requested Range Not Satisfiable")
            return None

        self.send_response(206)
        self.send_header("Content-type", ctype)

        if last is None or last >= file_len:
            last = file_len - 1
        response_length = last - first + 1

        self.send_header("Content-Range", "bytes %s-%s/%s" % (first, last, file_len))
        self.send_header("Content-Length", str(response_length))
        self.send_header("Last-Modified", self.date_time_string(fs.st_mtime))
        self.end_headers()
        return f

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        return SimpleHTTPRequestHandler.end_headers(self)

    def copyfile(self, source, outputfile):
        if not self.range:
            return SimpleHTTPRequestHandler.copyfile(self, source, outputfile)

        # SimpleHTTPRequestHandler uses shutil.copyfileobj, which doesn't let
        # you stop the copying before the end of the file.
        start, stop = self.range  # set in send_head()
        copy_byte_range(source, outputfile, start, stop)


class ThreadedTCPServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    allow_reuse_address = True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple Python Web Server with Range Support"
    )
    parser.add_argument(
        "--root",
        default=os.getcwd(),
        help="Root path to serve files from (default: current working directory)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to listen on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--bind", default="0.0.0.0", help="IP address to bind to (default: 0.0.0.0)"
    )
    args = parser.parse_args()

    handler = functools.partial(RangeRequestHandler, directory=args.root)

    webbrowser.open(f"http://{args.bind}:{args.port}")

    with ThreadedTCPServer((args.bind, args.port), handler) as httpd:
        print(
            f"Serving HTTP on {args.bind} port {args.port} (http://{args.bind}:{args.port}/)"
        )
        httpd.serve_forever()
