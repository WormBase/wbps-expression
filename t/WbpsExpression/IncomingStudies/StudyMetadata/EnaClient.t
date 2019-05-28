use strict;
use warnings;

use File::Temp qw/tempdir/;

use Test::More;
use Test::Mock::LWP;
use WbpsExpression::IncomingStudies::StudyMetadata::EnaClient;


my $study_id = "ERP006623";
my $bioproject_id = "PRJEB6948";
my $second_study_id  = "SRP124650";
my $second_bioproject_id = "PRJNA417697"; 

my ($study_xml, $bioproject_xml, $second_study_xml, $second_bioproject_xml) 
 = map {s/^\s*//; s/\s*$//; $_} split("---", do { local $/; <DATA> });

$Mock_ua->mock(get => sub { 
  our (undef, $url) = @_;

  $Mock_response->mock( decoded_content => sub { 
        if ( $url =~ /$study_id/ ) {
            return $study_xml;
        } elsif ($url =~ /$bioproject_id/){
            return $bioproject_xml;
        } elsif ( $url =~ /$second_study_id/ ) {
            return $second_study_xml;
        }
        elsif ($url =~ /$second_bioproject_id/ ) {
            return $second_bioproject_xml;
        }
        else {
            return undef;
        }
     } 
  );
  return $Mock_response;
});

my $root_dir          = tempdir( CLEANUP => 1 );
my $species           = "fasciola_hepatica";
my $assembly          = "Fasciola_10x_pilon";
my $rnaseqer_metadata = bless {
    metadata => {
        $assembly => {
            $study_id          => {},
            $second_study_id => {},
        }
    }
  },
  'PublicResources::Resources::RnaseqerMetadata';

is_deeply(WbpsExpression::IncomingStudies::StudyMetadata::EnaClient::get_study_metadata($study_id), 
  {
    'attributes' => {},
    'bioproject' => 'PRJEB6948',
    'description' => 'RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton\'s group and sequenced by Genome Quebec.',
    'ena_first_public' => '2014-12-31',
    'ena_last_update' => '2016-04-19',
    'pubmed_refs' => [
      '25887684',
      25887684
    ],
    'resource_links' => [
      [
        'ArrayExpress',
        'E-MTAB-451 in ArrayExpress',
        'http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-451'
      ]
    ],
    'submitting_centre' => 'University of Liverpool',
    'title' => 'Some RNA-seq reads from different developmental stages of the liver fluke Fasciola hepatica'
  }, "Test case from canned input");

is_deeply(WbpsExpression::IncomingStudies::StudyMetadata::EnaClient::get_study_metadata($second_study_id),
  {
    'attributes' => {},
    'bioproject' => 'PRJNA417697',
    'description' => '',
    'ena_first_public' => '2017-11-30',
    'ena_last_update' => '2017-11-30',
    'pubmed_refs' => [],
    'resource_links' => [],
    'submitting_centre' => 'Pharmacology, UT Southwestern',
    'title' => 'Transcriptional Profiling of Schistosoma mansoni f zfp-1-1(RNAi) parasites'
  }, "Second case from canned input");

done_testing();
__DATA__
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="ERP006623&amp;display=xml">
<STUDY alias="ena-STUDY-LIV-10-08-2014-19:12:12:707-88" center_name="University of Liverpool" accession="ERP006623">
     <IDENTIFIERS>
          <PRIMARY_ID>ERP006623</PRIMARY_ID>
          <SECONDARY_ID>PRJEB6948</SECONDARY_ID>
          <SUBMITTER_ID namespace="LIV">ena-STUDY-LIV-10-08-2014-19:12:12:707-88</SUBMITTER_ID>
     </IDENTIFIERS>
     <DESCRIPTOR>
          <STUDY_TITLE>Some RNA-seq reads form different developmental stages of the liver fluke Fasciola hepatica</STUDY_TITLE>
          <STUDY_ABSTRACT>RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.</STUDY_ABSTRACT>
          <STUDY_DESCRIPTION>RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.</STUDY_DESCRIPTION>
          <CENTER_PROJECT_NAME>Transcriptome of Fasciola hepatica</CENTER_PROJECT_NAME>
          <STUDY_TYPE existing_study_type="Other"/>
     </DESCRIPTOR>
     <STUDY_LINKS>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>PUBMED</DB>
                    <ID>25887684</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <URL_LINK>
                    <LABEL>E-MTAB-451 in ArrayExpress</LABEL>
                    <URL>http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-451</URL>
               </URL_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>ERS524684,ERS524692,ERS524694-ERS524695,ERS524697</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>ERX535559-ERX535563</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>ERR577156-ERR577160</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>ERA345947</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP006623&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERP006623&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
     </STUDY_LINKS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>232639591</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>46527918200</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2014-12-31</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2016-04-19</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>
---
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="PRJEB6948&amp;display=xml">
<PROJECT accession="PRJEB6948" alias="ena-STUDY-LIV-10-08-2014-19:12:12:707-88" center_name="University of Liverpool">
     <IDENTIFIERS>
          <PRIMARY_ID>PRJEB6948</PRIMARY_ID>
          <SECONDARY_ID>ERP006623</SECONDARY_ID>
          <SUBMITTER_ID namespace="LIV">ena-STUDY-LIV-10-08-2014-19:12:12:707-88</SUBMITTER_ID>
     </IDENTIFIERS>
     <NAME>Transcriptome of Fasciola hepatica</NAME>
     <TITLE>Some RNA-seq reads form different developmental stages of the liver fluke Fasciola hepatica</TITLE>
     <DESCRIPTION>RNA was prepared from various stages of the liver fluke Fasciola hepatica by John Dalton's group and sequenced by Genome Quebec.</DESCRIPTION>
     <SUBMISSION_PROJECT>
          <SEQUENCING_PROJECT/>
     </SUBMISSION_PROJECT>
     <PROJECT_LINKS>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>PUBMED</DB>
                    <ID>25887684</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>ERA345947</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>ERS524684,ERS524692,ERS524694-ERS524695,ERS524697</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>ERX535559-ERX535563</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>ERR577156-ERR577160</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB6948&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB6948&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </PROJECT_LINK>
     </PROJECT_LINKS>
     <PROJECT_ATTRIBUTES>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2014-12-31</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2016-05-20</VALUE>
          </PROJECT_ATTRIBUTE>
     </PROJECT_ATTRIBUTES>
</PROJECT>
</ROOT>
---
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="SRP124650&amp;display=xml">
<STUDY accession="SRP124650" alias="GSE106693" center_name="GEO" broker_name="NCBI">
     <IDENTIFIERS>
          <PRIMARY_ID>SRP124650</PRIMARY_ID>
          <SECONDARY_ID>PRJNA417697</SECONDARY_ID>
          <EXTERNAL_ID label="primary" namespace="BioProject">PRJNA417697</EXTERNAL_ID>
          <EXTERNAL_ID namespace="GEO">GSE106693</EXTERNAL_ID>
          <SUBMITTER_ID namespace="GEO">GSE106693</SUBMITTER_ID>
     </IDENTIFIERS>
     <DESCRIPTOR>
          <STUDY_TITLE>Transcriptional Profiling of Schistosoma mansoni f zfp-1-1(RNAi) parasites</STUDY_TITLE>
          <STUDY_TYPE existing_study_type="Transcriptome Analysis"/>
          <STUDY_ABSTRACT>We distrupted the expression of the Schistosoma mansoni zfp-1-1 gene using RNA interference and examined the transcriptional effects by illumina sequencing Overall design: The transcriptional effects of zfp-1-1 RNAi were examined by RNAseq</STUDY_ABSTRACT>
          <CENTER_PROJECT_NAME>GSE106693</CENTER_PROJECT_NAME>
     </DESCRIPTOR>
     <STUDY_LINKS>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>SRS2672175-SRS2672180</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>SRX3375885-SRX3375890</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>SRR6269807-SRR6269812</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>SRA629563</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP124650&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP124650&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
     </STUDY_LINKS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>175920283</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>13369941508</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2017-11-30</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2017-11-30</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>
---
<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="PRJNA417697&amp;display=xml">
<PROJECT accession="PRJNA417697" alias="PRJNA417697" center_name="Pharmacology, UT Southwestern">
     <IDENTIFIERS>
          <PRIMARY_ID>PRJNA417697</PRIMARY_ID>
          <SECONDARY_ID>SRP124650</SECONDARY_ID>
          <EXTERNAL_ID namespace="GEO">GSE106693</EXTERNAL_ID>
          <SUBMITTER_ID namespace="Pharmacology, UT Southwestern">PRJNA417697</SUBMITTER_ID>
     </IDENTIFIERS>
     <NAME>Transcriptional Profiling of Schistosoma mansoni f zfp-1-1(RNAi) parasites</NAME>
     <TITLE>Transcriptional Profiling of Schistosoma mansoni f zfp-1-1(RNAi) parasites</TITLE>
     <DESCRIPTION>We distrupted the expression of the Schistosoma mansoni zfp-1-1 gene using RNA interference and examined the transcriptional effects by illumina sequencing Overall design: The transcriptional effects of zfp-1-1 RNAi were examined by RNAseq</DESCRIPTION>
     <SUBMISSION_PROJECT>
          <SEQUENCING_PROJECT/>
          <ORGANISM>
               <TAXON_ID>6183</TAXON_ID>
               <SCIENTIFIC_NAME>Schistosoma mansoni</SCIENTIFIC_NAME>
          </ORGANISM>
     </SUBMISSION_PROJECT>
     <RELATED_PROJECTS>
          <RELATED_PROJECT>
               <PARENT_PROJECT accession="PRJNA12599"/>
          </RELATED_PROJECT>
     </RELATED_PROJECTS>
     <PROJECT_LINKS>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>SRS2672175-SRS2672180</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>SRX3375885-SRX3375890</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>SRR6269807-SRR6269812</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA417697&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA417697&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </PROJECT_LINK>
     </PROJECT_LINKS>
     <PROJECT_ATTRIBUTES>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-REFSEQ</TAG>
               <VALUE>N</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>PROJECT-ID</TAG>
               <VALUE>417697</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>NCBI-PROJECT-TYPE</TAG>
               <VALUE>SUBMISSION</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2017-11-10</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2018-01-28</VALUE>
          </PROJECT_ATTRIBUTE>
     </PROJECT_ATTRIBUTES>
</PROJECT>
</ROOT>
