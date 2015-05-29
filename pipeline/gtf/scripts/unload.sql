\set ON_ERROR_STOP on
\set ECHO all

BEGIN;

SELECT ge.id, ge.track_id INTO TEMPORARY TABLE genes_to_delete FROM
gene ge, coordinates co, chromosome ch, genome gg WHERE
ch.genome_id = gg.id AND
co.chromosome_id = ch.id AND
ge.track_id = co.id AND
gg.version = 'GRCh38.p2';

select count(*) from genes_to_delete;

SELECT id, track_id INTO TEMPORARY TABLE transcripts_to_delete FROM transcript WHERE
gene_id in (SELECT id from genes_to_delete);

SELECT id, track_id INTO TEMPORARY TABLE exons_to_delete FROM exon WHERE
transcript_id in (SELECT id from transcripts_to_delete);

SELECT id, track_id INTO TEMPORARY TABLE cds_to_delete FROM cds WHERE
gene_id in (SELECT id from genes_to_delete);

SELECT id, track_id INTO TEMPORARY TABLE cdsregions_to_delete FROM cdsregion WHERE
cds_id in (SELECT id from cds_to_delete);

DELETE FROM cdsregion WHERE id IN cdsregions_to_delete;
DELETE FROM cds WHERE id IN cds_to_delete;
DELETE FROM exon WHERE id IN exons_to_delete;
DELETE FROM transcript WHERE id IN transcripts_to_delete;
DELETE FROM gene WHERE id IN genes_to_delete;

DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM cdsregions_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM cds_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM exons_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM transcripts_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM genes_to_delete );

COMMIT;
