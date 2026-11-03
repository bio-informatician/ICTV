import importlib.util
import logging
from pathlib import Path

import vmr_ictv_merge
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "upsert_vjdb", ROOT / "m.upsert_vjdb.py"
)
upsert_vjdb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(upsert_vjdb)


def schema_fields():
    fields = {
        "insdc_accession": ("0aa", "keyword"),
        "is_refseq_sequence": ("0ai", "boolean"),
        "is_genbank_sequence": ("0ah", "boolean"),
        "sequence_source": ("0dk", "keyword"),
        "realm_ictv": ("0fi", "text"),
        "subrealm_ictv": ("0fj", "text"),
        "kingdom_ictv": ("0fk", "text"),
        "subkingdom_ictv": ("0fl", "text"),
        "phylum_ictv": ("0fm", "text"),
        "subphylum_ictv": ("0fn", "text"),
        "class_ictv": ("0fo", "text"),
        "subclass_ictv": ("0fp", "text"),
        "order_ictv": ("0fq", "text"),
        "suborder_ictv": ("0fr", "text"),
        "family_ictv": ("0fs", "text"),
        "subfamily_ictv": ("0ft", "text"),
        "genus_ictv": ("0fu", "text"),
        "subgenus_ictv": ("0fv", "text"),
        "species_ictv": ("0fw", "text"),
        "ictv_id": ("0gz", "keyword"),
        "is_ictv_exemplar": ("0gw", "boolean"),
        "host_group_ictv": ("0gy", "text"),
        "virus_name_ictv": ("0ha", "text"),
        "abbreviation": ("0cn", "keyword"),
        "isolate": ("0gn", "text"),
        "genome_completeness": ("0ao", "keyword"),
        "organism_name": ("0aq", "text"),
    }
    return {
        field_id: {"index": index, "type": field_type, "name": field_id}
        for field_id, (index, field_type) in fields.items()
    }


def existing_document():
    return {
        "0aa": ["AB123456.1"],
        "0dn": "vj000000000001",
        "0do": 1,
        "0cx": 1,
        "0hd": True,
        "0bg": "VMR-v1",
        "0fg": "2024-01-01",
        "0ai": False,
        "0ah": True,
        "0cd": False,
        "0hj": False,
        "0hl": False,
        "0dk": ["GenBank"],
    }


def merge_row(row, old_doc=None):
    schema = schema_fields()
    normalized = {}
    warnings = []
    vmr_ictv_merge.merge_vmr_into_normalized(
        normalized, [row], old_doc, schema, "accession AB123456", warnings
    )
    return normalized, warnings


class FakeResponse:
    def raise_for_status(self):
        pass

    def json(self):
        return {
            "hits": {
                "hits": [
                    {
                        "_id": "existing",
                        "_source": {
                            "0hd": True,
                            "0aa": ["ab123456.7"],
                            "0dn": "vj000000000001",
                        },
                    }
                ]
            }
        }


class FakeSession:
    def post(self, url, json, verify, timeout):
        self.request = (url, json, verify, timeout)
        return FakeResponse()


class SettingsResponse:
    def __init__(self, settings):
        self.settings = settings

    def raise_for_status(self):
        pass

    def json(self):
        return {"genome": {"settings": self.settings}}


class SettingsSession:
    def __init__(self, settings):
        self.settings = settings
        self.puts = []

    def get(self, url, params, verify, timeout):
        return SettingsResponse(self.settings)

    def put(self, url, json, verify, timeout):
        self.puts.append(json)
        return SettingsResponse({})


def test_bulk_errors_raise_to_fail_the_run():
    response = {
        "errors": True,
        "items": [{"index": {"error": {"type": "blocked", "reason": "read-only"}}}],
    }

    with pytest.raises(RuntimeError, match="metadata bulk failed"):
        upsert_vjdb.check_bulk_errors(response, logging.getLogger("test"), "metadata bulk")


def test_writable_genome_index_settings_are_not_changed(monkeypatch):
    session = SettingsSession({})
    expected = {"errors": False, "items": []}
    monkeypatch.setattr(upsert_vjdb, "execute_bulk", lambda *args, **kwargs: expected)

    result = upsert_vjdb.execute_genome_bulk(
        session, "https://example.invalid", "genome", [], True, 900_000
    )

    assert result is expected
    assert session.puts == []


def test_blocked_genome_index_is_restored_after_write(monkeypatch):
    session = SettingsSession({"index.blocks.write": "true"})
    monkeypatch.setattr(
        upsert_vjdb, "execute_bulk",
        lambda *args, **kwargs: {"errors": False, "items": []},
    )

    upsert_vjdb.execute_genome_bulk(
        session, "https://example.invalid", "genome", [], True, 900_000
    )

    assert session.puts == [
        {"index.blocks.write": None},
        {"index.blocks.write": True},
    ]


def test_accession_matching_uses_versionless_0aa_and_segmented_vmr_accessions():
    schema = schema_fields()
    session = FakeSession()
    db_lookup = upsert_vjdb.find_active_by_accession(
        session, "https://example.invalid", "metadata", schema,
        ["AB123456"], True,
    )
    assert db_lookup["AB123456"][0] == "existing"
    assert session.request[0] == "https://example.invalid/metadata/_search"
    assert session.request[1]["query"]["bool"]["filter"] == [
        {"term": {"0hd": True}}
    ]
    assert session.request[2:] == (True, 30)

    vmr_index, stats = vmr_ictv_merge.build_vmr_accession_index([
        {
            "Virus_GENBANK_accession":
                "DNA-A: ab123456.2; DNA-B: NC_123456.1",
            "Species": "Example species",
        }
    ])
    assert vmr_index["AB123456"][0]["Species"] == "Example species"
    assert "NC_123456" in vmr_index
    assert stats["distinct_base_accessions"] == 2


def test_vmr_columns_resolve_to_configured_schema_field_ids_and_indices():
    row = {
        "Realm": "Realm name",
        "Subrealm": "Subrealm name",
        "Kingdom": "Kingdom name",
        "Subkingdom": "Subkingdom name",
        "Phylum": "Phylum name",
        "Subphylum": "Subphylum name",
        "Class": "Class name",
        "Subclass": "Subclass name",
        "Order": "Order name",
        "Suborder": "Suborder name",
        "Family": "Family name",
        "Subfamily": "Subfamily name",
        "Genus": "Genus name",
        "Subgenus": "Subgenus name",
        "Species": "Species name",
        "ICTV_ID": "ICTV202600001",
        "Exemplar_or_additional_isolate": "E",
        "Host_source": "bacteria",
        "Virus_name_s_": "Example virus",
        "Virus_name_abbreviation_s_": "EXV",
        "Virus_isolate_designation": "isolate 1",
        "Genome_coverage": "Complete genome",
    }
    normalized, warnings = merge_row(row)
    assert not warnings
    assert normalized == {
        "realm_ictv": "Realm name",
        "subrealm_ictv": "Subrealm name",
        "kingdom_ictv": "Kingdom name",
        "subkingdom_ictv": "Subkingdom name",
        "phylum_ictv": "Phylum name",
        "subphylum_ictv": "Subphylum name",
        "class_ictv": "Class name",
        "subclass_ictv": "Subclass name",
        "order_ictv": "Order name",
        "suborder_ictv": "Suborder name",
        "family_ictv": "Family name",
        "subfamily_ictv": "Subfamily name",
        "genus_ictv": "Genus name",
        "subgenus_ictv": "Subgenus name",
        "species_ictv": "Species name",
        "ictv_id": "ICTV202600001",
        "is_ictv_exemplar": True,
        "host_group_ictv": "bacteria",
        "virus_name_ictv": "Example virus",
        "abbreviation": "EXV",
        "isolate": "isolate 1",
        "genome_completeness": "Complete genome",
    }

    errors = []
    mapped = upsert_vjdb.map_normalized_fields(
        normalized, schema_fields(), "test", errors, []
    )
    assert errors == []
    assert mapped["0fi"] == "Realm name"
    assert mapped["0fj"] == "Subrealm name"
    assert mapped["0fk"] == "Kingdom name"
    assert mapped["0fl"] == "Subkingdom name"
    assert mapped["0fm"] == "Phylum name"
    assert mapped["0fn"] == "Subphylum name"
    assert mapped["0fo"] == "Class name"
    assert mapped["0fp"] == "Subclass name"
    assert mapped["0fq"] == "Order name"
    assert mapped["0fr"] == "Suborder name"
    assert mapped["0fs"] == "Family name"
    assert mapped["0ft"] == "Subfamily name"
    assert mapped["0fu"] == "Genus name"
    assert mapped["0fv"] == "Subgenus name"
    assert mapped["0fw"] == "Species name"
    assert mapped["0gz"] == "ICTV202600001"
    assert mapped["0gw"] is True
    assert mapped["0gy"] == "bacteria"
    assert mapped["0ha"] == "Example virus"
    assert mapped["0cn"] == "EXV"
    assert mapped["0gn"] == "isolate 1"
    assert mapped["0ao"] == "Complete genome"


def test_existing_record_fills_empty_fields_and_preserves_blanks():
    row = {
        "Species": "New ICTV species",
        "Host_source": None,
        "Virus_name_s_": " ",
        "Genome_coverage": "",
    }
    old_doc = existing_document()
    old_doc.update({
        "0fw": "",
        "0gy": "existing host",
        "0ha": "existing virus name",
        "0ao": "Complete genome",
    })
    normalized, _ = merge_row(row, old_doc)

    result = upsert_vjdb.build_record_actions(
        "vj000000000001", "existing", old_doc, normalized, "AB123456.1",
        "", schema_fields(),
        {
            "version_change": {"method": "test", "reason": "test"},
            "content_update": {"method": "test", "reason": "test"},
        },
        "VMR-v1", "test",
    )

    assert result["event"] == "insertion"
    assert result["fill_update"]["doc"]["0fw"] == "New ICTV species"
    assert "0gy" not in result["fill_update"]["doc"]
    assert "0ha" not in result["fill_update"]["doc"]
    assert "0ao" not in result["fill_update"]["doc"]
    assert result["change_log_doc"] is None


def test_existing_record_changed_value_uses_content_update_and_change_log():
    old_doc = existing_document()
    old_doc["0gz"] = "ICTV202500001"
    normalized, _ = merge_row({"ICTV_ID": "ICTV202600001"}, old_doc)

    result = upsert_vjdb.build_record_actions(
        "vj000000000001", "existing", old_doc, normalized, "AB123456.1",
        "", schema_fields(),
        {
            "version_change": {"method": "test", "reason": "test"},
            "content_update": {"method": "test", "reason": "test"},
        },
        "VMR-v1", "test",
    )

    assert result["event"] == "content-update"
    assert result["new_meta_doc"]["0gz"] == "ICTV202600001"
    assert result["new_meta_doc"]["0do"] == 2
    assert result["change_log_doc"]["changes"] == [{
        "field": "0gz",
        "old_value": '"ICTV202500001"',
        "new_value": '"ICTV202600001"',
    }]


def test_version_diff_preserves_change_log_and_sequence_write_behavior():
    old_doc = existing_document()
    normalized, _ = merge_row({"ICTV_ID": "ICTV202600001"}, old_doc)

    result = upsert_vjdb.build_record_actions(
        "vj000000000001", "existing", old_doc, normalized, "AB123456.2",
        "ACGT", schema_fields(),
        {
            "version_change": {"method": "version test", "reason": "new accession version"},
            "content_update": {"method": "content test", "reason": "metadata changed"},
        },
        "VMR-v1", "test",
    )

    assert result["event"] == "version-diff"
    assert result["old_doc_update"] == {
        "_id": "existing",
        "doc": {"0hd": False},
    }
    assert result["new_meta_doc"]["0gz"] == "ICTV202600001"
    assert result["new_meta_doc"]["0do"] == 2
    assert result["new_meta_doc"]["0cx"] == 2
    assert result["change_log_doc"]["method"] == "version test"
    assert result["genome_doc"]["0aa"] == "AB123456.2"
    assert result["genome_doc"]["0cx"] == 2
    assert result["genome_doc"]["0hh"] == "ACGT"


def test_blank_vmr_value_does_not_clear_existing_field():
    old_doc = existing_document()
    old_doc["0gz"] = "ICTV202500001"
    normalized, _ = merge_row({"ICTV_ID": None}, old_doc)

    result = upsert_vjdb.build_record_actions(
        "vj000000000001", "existing", old_doc, normalized, "AB123456.1",
        "", schema_fields(),
        {
            "version_change": {"method": "test", "reason": "test"},
            "content_update": {"method": "test", "reason": "test"},
        },
        "VMR-v1", "test",
    )

    assert normalized == {}
    assert result["event"] == "insertion"
    assert "0gz" not in result["fill_update"]["doc"]
    assert old_doc["0gz"] == "ICTV202500001"


def test_new_record_contains_vmr_metadata_and_writes_fasta_sequence():
    normalized, _ = merge_row({
        "Realm": "Realm name",
        "Species": "Species name",
        "ICTV_ID": "ICTV202600001",
        "Exemplar_or_additional_isolate": "A",
    })
    normalized["virus_organism_name"] = "NCBI-resolved organism"

    result = upsert_vjdb.build_brand_new_insertion(
        "vj000000000002", normalized, "AB123456.1", "ACGT",
        schema_fields(), "VMR-v1", "new record",
    )

    assert result["errors"] == []
    assert result["new_meta_doc"]["0fi"] == "Realm name"
    assert result["new_meta_doc"]["0fw"] == "Species name"
    assert result["new_meta_doc"]["0gz"] == "ICTV202600001"
    assert result["new_meta_doc"]["0gw"] is False
    assert result["new_meta_doc"]["0aq"] == "NCBI-resolved organism"
    assert result["new_meta_doc"]["0dn"] == "vj000000000002"
    assert result["genome_doc"]["0aa"] == "AB123456.1"
    assert result["genome_doc"]["0hh"] == "ACGT"


def test_columns_without_database_destinations_are_reported():
    records = [{
        "Virus_GENBANK_accession": "AB123456",
        "Isolate_ID": "VMR123",
        "Species_Sort": 1,
        "Isolate_Sort": 1,
        "Genome": "dsDNA",
        "Accessions_Link": "NCBI Nucleotide",
        "from_url": "https://example.invalid/AB123456",
    }]
    columns, count = vmr_ictv_merge.report_unrecognized_keys(records)

    assert count == 6
    assert columns == [
        "Accessions_Link",
        "Genome",
        "Isolate_ID",
        "Isolate_Sort",
        "Species_Sort",
        "from_url",
    ]