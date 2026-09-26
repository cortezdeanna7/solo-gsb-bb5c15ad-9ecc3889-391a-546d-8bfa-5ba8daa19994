"""Tests for the /partial-digest endpoint and its exact probability logic."""

from fastapi.testclient import TestClient

from app.main import app
from app.partial import (
    PartialEnzyme,
    circular_expectations,
    effective_cut_positions,
    linear_expectations,
)

client = TestClient(app)


# ---------- core logic ----------


def test_effective_positions_merge_shared_cuts():
    # Two enzymes cutting the same two positions with p=0.5 each:
    # nick probability per position = 1 - (1-0.5)^2 = 0.75.
    enzymes = [
        PartialEnzyme(name="A", site="GAATTC", cut_offset=1, p=0.5),
        PartialEnzyme(name="B", site="GAATTC", cut_offset=1, p=0.5),
    ]
    positions, probs = effective_cut_positions("AGAATTCT", enzymes, circular=False)
    assert positions == [2, 6]
    assert probs == [0.75, 0.75]


def test_effective_positions_drop_zero_probability():
    enzymes = [PartialEnzyme(name="A", site="CG", cut_offset=1, p=0.0)]
    positions, probs = effective_cut_positions("ACGT", enzymes, circular=False)
    assert positions == [] and probs == []


def test_effective_positions_exclude_linear_ends():
    # k=0 puts both the forward and the reverse-complement cut on the
    # molecule ends (0 and 6); linear ends are not random nicks.
    enzymes = [PartialEnzyme(name="A", site="GAATTC", cut_offset=0, p=1.0)]
    positions, _ = effective_cut_positions("GAATTC", enzymes, circular=False)
    assert positions == []
    # On a circle the same cuts collapse (mod 6) to one real position.
    positions, probs = effective_cut_positions("GAATTC", enzymes, circular=True)
    assert positions == [0]
    assert probs == [1.0]


def test_linear_expectations_single_cut():
    expected, no_cut = linear_expectations(4, [2], [0.5])
    assert expected == {2: 1.0, 4: 0.5}  # same-length fragments accumulate
    assert no_cut == 0.5


def test_linear_expectations_certain_cuts():
    # p=1 everywhere: only the deterministic fragments remain.
    expected, no_cut = linear_expectations(20, [5, 9, 13, 17], [1.0] * 4)
    assert expected == {3: 1.0, 4: 3.0, 5: 1.0}
    assert no_cut == 0.0


def test_circular_expectations_single_cut_full_length():
    # Zero or exactly one nick both yield one full-length fragment.
    expected, no_cut = circular_expectations(6, [0], [0.5])
    assert expected == {6: 1.0}
    assert no_cut == 0.5  # only the zero-nick case counts


def test_circular_expectations_two_cuts():
    expected, no_cut = circular_expectations(8, [2, 6], [0.5, 0.5])
    # arcs of 4 with prob 0.25 each way; full length with 0.25 + 0.5.
    assert expected == {4: 0.5, 8: 0.75}
    assert no_cut == 0.25


def test_circular_expectations_certain_cuts():
    expected, no_cut = circular_expectations(6, [1, 5], [1.0, 1.0])
    assert expected == {2: 1.0, 4: 1.0}
    assert no_cut == 0.0


# ---------- HTTP API ----------


def test_partial_digest_linear_certain_matches_deterministic():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGTGAATTCGGGAATTCAA",
            "topology": "linear",
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [
            {"length": 3, "expected_count": 1.0},
            {"length": 4, "expected_count": 3.0},
            {"length": 5, "expected_count": 1.0},
        ],
        "no_cut_probability": 0.0,
    }


def test_partial_digest_linear_four_positions():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGTGAATTCGGGAATTCAA",
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 0.5}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [
            {"length": 3, "expected_count": 0.5},
            {"length": 4, "expected_count": 0.75},
            {"length": 5, "expected_count": 0.5},
            {"length": 7, "expected_count": 0.25},
            {"length": 8, "expected_count": 0.25},
            {"length": 9, "expected_count": 0.25},
            {"length": 11, "expected_count": 0.125},
            {"length": 12, "expected_count": 0.0625},
            {"length": 13, "expected_count": 0.125},
            {"length": 15, "expected_count": 0.0625},
            {"length": 17, "expected_count": 0.0625},
            {"length": 20, "expected_count": 0.0625},
        ],
        "no_cut_probability": 0.0625,
    }


def test_partial_digest_shared_position_two_enzymes():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "AGAATTCT",
            "enzymes": [
                {"name": "EcoRI", "site": "GAATTC", "cut_offset": 1, "p": 0.5},
                {"name": "EcoRI-like", "site": "GAATTC", "cut_offset": 1, "p": 0.5},
            ],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [
            {"length": 2, "expected_count": 1.5},
            {"length": 4, "expected_count": 0.5625},
            {"length": 6, "expected_count": 0.375},
            {"length": 8, "expected_count": 0.0625},
        ],
        "no_cut_probability": 0.0625,
    }


def test_partial_digest_circular_single_cut():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTC",
            "topology": "circular",
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 0, "p": 0.5}],
        },
    )
    assert resp.status_code == 200
    # Zero nicks (0.5) and exactly one nick (0.5) both give one full-length
    # fragment; only the zero-nick case counts towards no_cut_probability.
    assert resp.json() == {
        "fragments": [{"length": 6, "expected_count": 1.0}],
        "no_cut_probability": 0.5,
    }


def test_partial_digest_circular_two_positions():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGTACGT",
            "topology": "circular",
            "enzymes": [{"name": "CgrI", "site": "CG", "cut_offset": 1, "p": 0.5}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [
            {"length": 4, "expected_count": 0.5},
            {"length": 8, "expected_count": 0.75},
        ],
        "no_cut_probability": 0.25,
    }


def test_partial_digest_empty_enzyme_list_means_no_cut():
    resp = client.post("/partial-digest", json={"sequence": "ACGT", "enzymes": []})
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [{"length": 4, "expected_count": 1.0}],
        "no_cut_probability": 1.0,
    }


def test_partial_digest_zero_probability_means_no_cut():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "CgrI", "site": "CG", "cut_offset": 1, "p": 0}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [{"length": 4, "expected_count": 1.0}],
        "no_cut_probability": 1.0,
    }


def test_partial_digest_linear_ends_are_not_random_nicks():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "GAATTC",
            "topology": "linear",
            "enzymes": [{"name": "EcoRI", "site": "GAATTC", "cut_offset": 0, "p": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [{"length": 6, "expected_count": 1.0}],
        "no_cut_probability": 1.0,
    }


def test_partial_digest_defaults_to_linear_and_lowercase_ok():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "acgt",
            "enzymes": [{"name": "CgrI", "site": "cg", "cut_offset": 1, "p": 0.5}],
        },
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [
            {"length": 2, "expected_count": 1.0},
            {"length": 4, "expected_count": 0.5},
        ],
        "no_cut_probability": 0.5,
    }


def test_partial_digest_max_length_sequence_ok():
    resp = client.post(
        "/partial-digest", json={"sequence": "A" * 2000, "enzymes": []}
    )
    assert resp.status_code == 200
    assert resp.json() == {
        "fragments": [{"length": 2000, "expected_count": 1.0}],
        "no_cut_probability": 1.0,
    }


# ---------- validation (HTTP 422) ----------


def test_reject_probability_out_of_range():
    for bad_p in (-0.1, 1.5):
        resp = client.post(
            "/partial-digest",
            json={
                "sequence": "ACGT",
                "enzymes": [{"name": "X", "site": "CG", "cut_offset": 1, "p": bad_p}],
            },
        )
        assert resp.status_code == 422


def test_reject_missing_probability():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "CG", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 422


def test_reject_sequence_longer_than_2000bp():
    resp = client.post(
        "/partial-digest", json={"sequence": "A" * 2001, "enzymes": []}
    )
    assert resp.status_code == 422


def test_reject_duplicate_enzyme_names():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [
                {"name": "X", "site": "AC", "cut_offset": 0, "p": 0.5},
                {"name": "X", "site": "GT", "cut_offset": 1, "p": 0.5},
            ],
        },
    )
    assert resp.status_code == 422


def test_reject_non_acgt_sequence_and_non_iupac_site():
    enzyme = {"name": "X", "site": "AC", "cut_offset": 0, "p": 0.5}
    resp = client.post(
        "/partial-digest", json={"sequence": "ACGTX", "enzymes": [enzyme]}
    )
    assert resp.status_code == 422
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "ACZ", "cut_offset": 0, "p": 0.5}],
        },
    )
    assert resp.status_code == 422


def test_reject_bad_offset_and_site_longer_than_sequence():
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "ACG", "cut_offset": 4, "p": 0.5}],
        },
    )
    assert resp.status_code == 422
    resp = client.post(
        "/partial-digest",
        json={
            "sequence": "ACGT",
            "enzymes": [{"name": "X", "site": "GAATTC", "cut_offset": 1, "p": 0.5}],
        },
    )
    assert resp.status_code == 422


# ---------- /digest and /design-panel regression ----------


def test_digest_endpoint_still_works_without_p():
    resp = client.post(
        "/digest",
        json={
            "sequence": "AGGATCCGAATTCA",
            "enzymes": [{"name": "BamHI", "site": "GGATCC", "cut_offset": 1}],
        },
    )
    assert resp.status_code == 200
    assert resp.json()["cuts"] == [{"name": "BamHI", "cuts": [2, 6]}]
