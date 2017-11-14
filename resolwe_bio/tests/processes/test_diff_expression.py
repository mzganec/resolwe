# pylint: disable=missing-docstring
from resolwe.flow.models import Data, Process
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class DiffExpProcessorTestCase(BioProcessTestCase):

    @tag_process('cuffdiff')
    def test_cuffdiff(self):
        with self.preparation_stage():
            inputs = {'src': 'cuffquant_1.cxb', 'source': 'UCSC'}
            cuffquant = self.run_process("upload-cxb", inputs)

            inputs = {'src': 'cuffquant_2.cxb', 'source': 'UCSC'}
            cuffquant2 = self.run_process("upload-cxb", inputs)

            annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz', source='UCSC')

        inputs = {
            'case': [cuffquant.id],
            'control': [cuffquant2.id],
            'annotation': annotation.id}
        cuffdiff = self.run_process('cuffdiff', inputs)
        self.assertFile(cuffdiff, 'raw', 'raw_cuffdiff.tab.gz', compression='gzip')
        self.assertFile(cuffdiff, 'de_file', 'de_file_cuffdiff.tab.gz', compression='gzip')
        self.assertJSON(cuffdiff, cuffdiff.output['de_json'], '', 'cuffdiff.json.gz')

    @tag_process('differentialexpression-bcm')
    def test_bayseq_bcm(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

            mappa = self.run_process("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'name': "00vs20",
            'case': [expression_1.id],
            'control': [expression_2.id],
            'replicates': ['1', '2'],
            'mappability': mappa.id}
        diff_exp = self.run_process('differentialexpression-bcm', inputs)
        self.assertJSON(diff_exp, diff_exp.output['volcano_plot'], '', 'bayseq_volcano.json.gz')

    @tag_process('differentialexpression-deseq2')
    def test_deseq2_genes(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', source='DICTYBASE')
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', source='DICTYBASE')
            expression_3 = self.prepare_expression(f_rc='exp_3_rc.tab.gz', source='DICTYBASE')
            expression_4 = self.prepare_expression(f_rc='exp_4_rc.tab.gz', source='DICTYBASE')

        inputs = {
            'case': [expression_1.pk, expression_3.pk],
            'control': [expression_2.pk, expression_4.pk],
            'filter': 0,
        }

        diff_exp = self.run_process('differentialexpression-deseq2', inputs)

        self.assertFileExists(diff_exp, 'raw')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'deseq2.json.gz')
        self.assertFields(diff_exp, 'source', 'DICTYBASE')

    @tag_process('differentialexpression-deseq2')
    def test_deseq2_transcripts(self):
        with self.preparation_stage():
            # Mock RSEM process
            process = Process.objects.create(
                name='RSEM mock process',
                requirements={
                    'expression-engine': 'jinja',
                    'resources': {
                        'network': True,
                    },
                    'executor': {
                        'docker': {
                            'image': 'resolwebio/base:ubuntu-17.04',
                        },
                    },
                },
                contributor=self.contributor,
                type='data:expression:rsem:',
                input_schema=[
                    {
                        'name': 'transcripts',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'transcripts',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'source',
                        'type': 'basic:string:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': r"""
re-import {{ transcripts.file_temp|default(transcripts.file) }} {{ transcripts.file }} "tab|gz" "tab" 1.0 compress
re-save-file transcripts "${NAME}.tab.gz"
re-save source DICTYBASE
"""
                }
            )

            inputs1 = {'transcripts': 'rsem_transcripts1.tab.gz'}
            rsem1 = self.run_process(process.slug, inputs1)
            self.assertFile(rsem1, 'transcripts', 'rsem_transcripts1.tab.gz', compression='gzip')

            inputs2 = {'transcripts': 'rsem_transcripts2.tab.gz'}
            rsem2 = self.run_process(process.slug, inputs2)
            self.assertFile(rsem2, 'transcripts', 'rsem_transcripts2.tab.gz', compression='gzip')

            inputs3 = {'transcripts': 'rsem_transcripts3.tab.gz'}
            rsem3 = self.run_process(process.slug, inputs3)
            self.assertFile(rsem3, 'transcripts', 'rsem_transcripts3.tab.gz', compression='gzip')

            inputs4 = {'transcripts': 'rsem_transcripts4.tab.gz'}
            rsem4 = self.run_process(process.slug, inputs4)
            self.assertFile(rsem4, 'transcripts', 'rsem_transcripts4.tab.gz', compression='gzip')

        inputs = {
            'case': [
                rsem1.pk,
                rsem4.pk,
            ],
            'control': [
                rsem2.pk,
                rsem3.pk,
            ],
            'filter': 0,
        }
        de_rsem = self.run_process('differentialexpression-deseq2', inputs)
        self.assertFile(de_rsem, 'raw', 'deseq_rsem.tab.gz', compression='gzip')
        self.assertJSON(de_rsem, de_rsem.output['de_json'], '', 'deseq_rsem.json.gz')
        self.assertFields(de_rsem, 'source', 'DICTYBASE')

    @tag_process('differentialexpression-deseq2')
    def test_deseq2_source(self):
        with self.preparation_stage():
            expression_dictybase = self.prepare_expression(source='DICTYBASE')
            expression_ucsc = self.prepare_expression(source='UCSC')

        inputs = {
            'case': [expression_dictybase.pk],
            'control': [expression_ucsc.pk]
        }

        self.run_process('differentialexpression-deseq2', inputs, Data.STATUS_ERROR)

    @tag_process('differentialexpression-limma')
    def test_limma(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_exp='exp_limma_1.tab.gz', f_type="Log2")
            expression_2 = self.prepare_expression(f_exp='exp_limma_2.tab.gz', f_type="Log2")
            expression_3 = self.prepare_expression(f_exp='exp_limma_3.tab.gz', f_type="Log2")
            expression_4 = self.prepare_expression(f_exp='exp_limma_4.tab.gz', f_type="Log2")

        inputs = {
            'case': [expression_1.pk, expression_2.pk],
            'control': [expression_3.pk, expression_4.pk]
        }

        diff_exp = self.run_process('differentialexpression-limma', inputs)
        self.assertFile(diff_exp, "raw", 'diffexp_limma.tab.gz', compression='gzip')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'limma.json.gz')

    @tag_process('differentialexpression-edger')
    def test_edger(self):
        with self.preparation_stage():
            inputs = {'rc': 'exp_1_rc.tab.gz', 'exp_name': 'Expression', 'source': 'DICTYBASE'}
            expression_1 = self.run_process('upload-expression', inputs)

            inputs = {'rc': 'exp_2_rc.tab.gz', 'exp_name': 'Expression', 'source': 'DICTYBASE'}
            expression_2 = self.run_process('upload-expression', inputs)

            inputs = {'rc': 'exp_3_rc.tab.gz', 'exp_name': 'Expression', 'source': 'DICTYBASE'}
            expression_3 = self.run_process('upload-expression', inputs)

            inputs = {'rc': 'exp_4_rc.tab.gz', 'exp_name': 'Expression', 'source': 'DICTYBASE'}
            expression_4 = self.run_process('upload-expression', inputs)

        inputs = {
            'case': [expression_1.pk, expression_3.pk],
            'control': [expression_2.pk, expression_4.pk]
        }

        diff_exp = self.run_process('differentialexpression-edger', inputs)
        self.assertFile(diff_exp, 'raw', 'diffexp_edgeR.tab.gz', compression='gzip')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'edgeR.json.gz')
        self.assertFields(diff_exp, 'source', 'DICTYBASE')
