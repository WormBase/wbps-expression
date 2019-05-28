use strict;
use warnings;
use Test::More;
use Test::Mock::LWP;
use WbpsExpression::IncomingStudies;
use File::Temp qw/tempdir/;
use File::Path qw/remove_tree/;


my ($runs_by_organism, $ftp_path_content, $study_metadata, $bioproject_metadata, $geo_no_results, $run_metadata) = map {s/^\s*//; s/\s*$//; $_} split("---", do { local $/; <DATA> });

$Mock_ua->mock(get => sub {
  our (undef, $url) = @_;
  $Mock_response->mock( decoded_content => sub {
    return $runs_by_organism if $url =~ m{getRunsByOrganism};
    return $ftp_path_content if $url =~ m{ftp://.*data/atlas/rnaseq.*/};
    return $study_metadata   if $url =~ m{ena/data/view/SRP070744};
    return $bioproject_metadata   if $url =~ m{ena/data/view/PRJNA312925};
    return $geo_no_results if $url =~ m{eutils};
    return $run_metadata if $url =~ m{getSampleAttributesPerRunByStudy};
    die $url;
  });
  return $Mock_response;
});


my $tmp = tempdir(CLEANUP => 1);
my $species = "anisakis_simplex";
my $assembly = "A_simplex_v1_5_4";
is_deeply(
  [[],[]],
  [WbpsExpression::IncomingStudies::update_studies($tmp, $species, "Different assembly")],
  "No data for different assembly"
);
my ($studies, $other_studies) = WbpsExpression::IncomingStudies::update_studies($tmp, $species, $assembly);
is_deeply(
  [[],["SRP070744"]],
  [$studies, [map {$_->{study_id}} @{$other_studies//[]}]],
  "Reject the study",
) or diag explain $studies, $other_studies;

remove_tree $tmp;

$WbpsExpression::IncomingStudies::exceptions{SRP070744} = "Few runs, still a great test case";
($studies, $other_studies) = WbpsExpression::IncomingStudies::update_studies($tmp, $species, $assembly);
is_deeply(
  [["SRP070744"], []],
  [[map {$_->{study_id}} @{$studies//[]}], $other_studies],
  "Upon adding an exception, decide this is after all a good study",
) or diag explain $studies, $other_studies;

done_testing;

__DATA__
[{"STUDY_ID":"SRP070744","SAMPLE_IDS":"SAMN04510371","BIOREP_ID":"SRR3184913","RUN_IDS":"SRR3184913","ORGANISM":"anisakis_simplex","REFERENCE_ORGANISM":"anisakis_simplex","STATUS":"Complete","ASSEMBLY_USED":"A_simplex_v1_5_4","ENA_LAST_UPDATED":"Tue Jun 28 2016 01:08:48","LAST_PROCESSED_DATE":"Fri Feb 09 2018 02:21:59","CRAM_LOCATION":"ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR318/003/SRR3184913/SRR3184913.cram","BEDGRAPH_LOCATION":"ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR318/003/SRR3184913/SRR3184913.bedgraph","BIGWIG_LOCATION":"ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/SRR318/003/SRR3184913/SRR3184913.bw","MAPPING_QUALITY":75}]
---
irap.versions.tsv
qc
SRR3184913.bedgraph
SRR3184913.bw
SRR3184913.cram
SRR3184913.cram.crai
SRR3184913.cram.md5
SRR3184913.nospliced.bedgraph
SRR3184913.nospliced.bw
SRR3184913.pe.exons.fpkm.dexseq.irap.tsv
SRR3184913.pe.exons.raw.dexseq.tsv
SRR3184913.pe.exons.tpm.dexseq.irap.tsv
SRR3184913.pe.genes.fpkm.htseq2.irap.tsv
SRR3184913.pe.genes.fpkm.kallisto.irap.tsv
SRR3184913.pe.genes.raw.htseq2.tsv
SRR3184913.pe.genes.raw.kallisto.tsv
SRR3184913.pe.genes.tpm.htseq2.irap.tsv
SRR3184913.pe.genes.tpm.kallisto.irap.tsv
SRR3184913.pe.hits.bam.gene.stats
SRR3184913.pe.hits.bam.stats
SRR3184913.pe.hits.bam.stats.csv
SRR3184913.pe.transcripts.dt.kallisto.irap.tsv
SRR3184913.pe.transcripts.fpkm.kallisto.irap.tsv
SRR3184913.pe.transcripts.raw.kallisto.tsv
SRR3184913.pe.transcripts.riu.kallisto.irap.tsv
SRR3184913.pe.transcripts.tpm.kallisto.irap.tsv
---

<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="SRP070744&amp;display=xml">
<STUDY accession="SRP070744" alias="PRJNA312925" center_name="BioProject" broker_name="NCBI">
     <IDENTIFIERS>
          <PRIMARY_ID>SRP070744</PRIMARY_ID>
          <SECONDARY_ID>PRJNA312925</SECONDARY_ID>
          <EXTERNAL_ID label="primary" namespace="BioProject">PRJNA312925</EXTERNAL_ID>
     </IDENTIFIERS>
     <DESCRIPTOR>
          <STUDY_TITLE>Detection and characterisation of putative allergens in Anisakis food-borne parasites using advanced transcriptomic and bioinformatic technologies</STUDY_TITLE>
          <STUDY_TYPE existing_study_type="Other"/>
          <STUDY_ABSTRACT>Background: Food-borne nematodes of the genus Anisakis are responsible for a widerange of illnesses (= anisakiasis), from self-limiting gastrointestinal forms to severesystemic allergic reactions, which are often misdiagnosed and under-reported. In orderto enhance and refine current diagnostic tools for anisakiasis, knowledge of the wholespectrum of parasite molecules acting as potential allergens is necessary.Methodology/Principal Findings: In this study, we employ high-throughput (Illumina)sequencing and bioinformatics technologies to characterise the transcriptomes of twoAnisakis species, A. simplex and A. pegreffii, and mine these annotated datasets tocompile lists of potential allergens from these parasites. A total of ~65,000,000 readswere generated from cDNA libraries for each species, and assembled into ~34,000transcripts (= Unigenes); ~18,000 peptides were predicted from each cDNA library andclassified based on homology searches, protein motifs and gene ontology andbiological pathway mapping. Using comparative analyses with sequence data availablein public databases, 36 (A. simplex) and 29 (A. pegreffii) putative allergens wereidentified, including sequences encoding 'novel' Anisakis allergenic proteins (i.e.cyclophilins and ABA-1 domain containing proteins).Conclusions/Significance: This study represents a first step towards providing theresearch community with a curated dataset to use as a molecular resource for futureinvestigations of poorly known putative Anisakis allergens, using functional genomics,proteomics and immunological tools. Ultimately, an improved knowledge of thebiological functions of these molecules in the parasite, as well as of their immunogenicproperties, will assist the development of comprehensive, reliable and robustdiagnostic tools.</STUDY_ABSTRACT>
          <CENTER_PROJECT_NAME>Anisakis spp.</CENTER_PROJECT_NAME>
     </DESCRIPTOR>
     <STUDY_LINKS>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>SRS1308153,SRS1308156</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>SRX1598161,SRX1598166</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>SRR3184913-SRR3184914</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>SRA357242</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP070744&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP070744&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
     </STUDY_LINKS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>91875463</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>18375092600</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2016-06-28</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2016-06-28</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>
---

<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="PRJNA312925&amp;display=xml">
<PROJECT accession="PRJNA312925" alias="PRJNA312925" center_name="University of Cambridge">
     <IDENTIFIERS>
          <PRIMARY_ID>PRJNA312925</PRIMARY_ID>
          <SECONDARY_ID>SRP070744</SECONDARY_ID>
          <SUBMITTER_ID namespace="University of Cambridge">PRJNA312925</SUBMITTER_ID>
     </IDENTIFIERS>
     <NAME>Anisakis spp.</NAME>
     <TITLE>Detection and characterisation of putative allergens in Anisakis food-borne parasites using advanced transcriptomic and bioinformatic technologies</TITLE>
     <DESCRIPTION><![CDATA[Background: Food-borne nematodes of the genus Anisakis are responsible for a wide</p><p>range of illnesses (= anisakiasis), from self-limiting gastrointestinal forms to severe</p><p>systemic allergic reactions, which are often misdiagnosed and under-reported. In order</p><p>to enhance and refine current diagnostic tools for anisakiasis, knowledge of the whole</p><p>spectrum of parasite molecules acting as potential allergens is necessary.</p><p>Methodology/Principal Findings: In this study, we employ high-throughput (Illumina)</p><p>sequencing and bioinformatics technologies to characterise the transcriptomes of two</p><p>Anisakis species, A. simplex and A. pegreffii, and mine these annotated datasets to</p><p>compile lists of potential allergens from these parasites. A total of ~65,000,000 reads</p><p>were generated from cDNA libraries for each species, and assembled into ~34,000</p><p>transcripts (= Unigenes); ~18,000 peptides were predicted from each cDNA library and</p><p>classified based on homology searches, protein motifs and gene ontology and</p><p>biological pathway mapping. Using comparative analyses with sequence data available</p><p>in public databases, 36 (A. simplex) and 29 (A. pegreffii) putative allergens were</p><p>identified, including sequences encoding 'novel' Anisakis allergenic proteins (i.e.</p><p>cyclophilins and ABA-1 domain containing proteins).</p><p>Conclusions/Significance: This study represents a first step towards providing the</p><p>research community with a curated dataset to use as a molecular resource for future</p><p>investigations of poorly known putative Anisakis allergens, using functional genomics,</p><p>proteomics and immunological tools. Ultimately, an improved knowledge of the</p><p>biological functions of these molecules in the parasite, as well as of their immunogenic</p><p>properties, will assist the development of comprehensive, reliable and robust</p><p>diagnostic tools.]]></DESCRIPTION>
     <SUBMISSION_PROJECT>
          <SEQUENCING_PROJECT/>
     </SUBMISSION_PROJECT>
     <PROJECT_LINKS>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>SRS1308153,SRS1308156</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>SRX1598161,SRX1598166</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>SRR3184913-SRR3184914</ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA312925&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </PROJECT_LINK>
          <PROJECT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJNA312925&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
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
               <VALUE>312925</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>NCBI-PROJECT-TYPE</TAG>
               <VALUE>SUBMISSION</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2016-06-29</VALUE>
          </PROJECT_ATTRIBUTE>
          <PROJECT_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2018-01-28</VALUE>
          </PROJECT_ATTRIBUTE>
     </PROJECT_ATTRIBUTES>
</PROJECT>
</ROOT>
---
<eSearchResult><Count>0</Count><RetMax>0</RetMax><RetStart>0</RetStart><QueryKey>1</QueryKey><WebEnv>NCID_1_167551032_130.14.22.76_9001_1558634035_1430741919_0MetA0_S_MegaStore</WebEnv><IdList/><TranslationSet/><QueryTranslation>(SRP070744[accn])</QueryTranslation><ErrorList><PhraseNotFound>SRP070744[accn]</PhraseNotFound></ErrorList><WarningList><OutputMessage>No items found.</OutputMessage></WarningList></eSearchResult>
---
[{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"Sample Name","VALUE":"Transcriptome of Anisakis simplex third stage larva","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"biomaterial provider","VALUE":"Hiromu Sugiyama","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"development stage","VALUE":"Third stage larva","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"geographic location","VALUE":"Japan: Tokyo","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"host","VALUE":"Scomber japonicus Scomber australasicus Trachurus japonicus","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"identified by","VALUE":"Fiona J. Baird","EFO_URL":"NA"},{"STUDY_ID":"SRP070744","RUN_ID":"SRR3184913","TYPE":"organism part","VALUE":"Whole larva","EFO_URL":"NA"}]
