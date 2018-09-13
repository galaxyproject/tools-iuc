class Wiggle:

    def fixedStepParser(self, line):
        value = line.strip()
        start_position = self.stepIdx * self.parserConfig['step'] + self.parserConfig['start']
        stop_position = start_position + self.parserConfig['span'] - 1
        self.stepIdx += 1

        for position in range(start_position, stop_position):
            yield (self.parserConfig['chrom'], position, value)

    def variableStepParser(self, line):
        (start, value) = line.strip().split()
        start = int(start)
        start_position = start
        stop_position = start + self.parserConfig['span']

        for position in range(start_position, stop_position):
            yield (self.parserConfig['chrom'], position, value)

    def walk(self, handle):

        parser = None
        for line in handle:
            if line.startswith('track'):
                continue
            elif line.startswith('fixedStep'):
                parser = self.fixedStepParser
                lineData = line.split()
                fields = {x.split('=')[0]: x.split('=')[1] for x in lineData[1:]}
                self.parserConfig = fields

                for numField in ('step', 'start', 'span'):
                    if numField in self.parserConfig:
                        self.parserConfig[numField] = int(self.parserConfig[numField])
                self.stepIdx = 0
            elif line.startswith('variableStep'):
                parser = self.variableStepParser
                lineData = line.split()
                fields = {x.split('=')[0]: x.split('=')[1] for x in lineData[1:]}
                # Default value
                if 'span' not in fields:
                    fields['span'] = 1
                self.parserConfig = fields

                for numField in ('span',):
                    if numField in self.parserConfig:
                        self.parserConfig[numField] = int(self.parserConfig[numField])

                self.stepIdx = 0
            elif len(line.strip()) == 0:
                continue
            else:
                for data in parser(line):
                    yield data
