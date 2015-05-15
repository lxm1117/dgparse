BEGIN;

SELECT track_id INTO TEMPORARY TABLE cdsregion_track_ids_to_delete FROM cdsregion;
SELECT track_id INTO TEMPORARY TABLE cds_track_ids_to_delete FROM cds;
SELECT track_id INTO TEMPORARY TABLE exon_track_ids_to_delete FROM exon;
SELECT track_id INTO TEMPORARY TABLE transcript_track_ids_to_delete FROM transcript;
SELECT track_id INTO TEMPORARY TABLE gene_track_ids_to_delete FROM gene;

DELETE FROM cdsregion;
DELETE FROM cds;
DELETE FROM exon;
DELETE FROM transcript;
DELETE FROM gene;

DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM exon_track_ids_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM transcript_track_ids_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM cdsregion_track_ids_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM cds_track_ids_to_delete );
DELETE FROM coordinates WHERE id IN ( SELECT track_id FROM gene_track_ids_to_delete );

COMMIT;

VACUUM ANALYZE;
