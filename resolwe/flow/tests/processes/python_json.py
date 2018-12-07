from resolwe.process import *


class PythonProcess(Process):
    """This is a process description."""
    slug = 'test-python-process-json'
    name = "Test Python Process that uses JSNO field"
    process_type = 'data:python'
    version = '0.1.2'
    requirements = {
        'executor': {
            'docker': {
                'image': 'resolwe/base:ubuntu-18.04',
            }
        }
    }

    class Input:
        """Input fields."""
        data = DataField('test', label="My input data")

    def run(self, inputs, outputs):
        print('Input data:', inputs.data)
        print('Input data:', inputs.data.storage)
