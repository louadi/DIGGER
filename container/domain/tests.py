import pytest
# Set up according to https://djangostars.com/blog/django-pytest-testing/ and
# https://stackoverflow.com/questions/60369047/pytest-using-parametized-fixture-vs-pytest-mark-parametrize

from .Process import nease_output as no

# Run with the command pytest, add -s if you want prints (might not work now with pytest-xdist), add -n auto if you want to run in parallel


input_files = [("Standard", "test_file1.json"), ("MAJIQ", "test_file2.json"), ("Whippet", "test_file2.json"),
               ("rMATS", "test_file2.json")]
organism = ["human", "mouse"]
predicted_ddis = [[], ["high"], ["mid"], ["low"], ["high", "mid"], ["mid", "low"],  ["high", "low"], ["high", "mid", "low"]]
p_values = [10, 1, 0.05, 0, -1]
deltas = [10, 1, 0.05, 0, -1]
majiq_conf_values = [10, 1, 0.05, 0, -1]
only_ddis_values = [True, False]
remove_not_in_frame_values = [True, False]
only_divisible_by_three_values = [True, False]


# Use pytest.mark.parametrize to generate all combinations of the parameters
@pytest.mark.parametrize("input", input_files)
@pytest.mark.parametrize("org", organism)
@pytest.mark.parametrize("pred_ddi", predicted_ddis)
@pytest.mark.parametrize("p_value", p_values)
@pytest.mark.parametrize("delta", deltas)
@pytest.mark.parametrize("majiq_conf", majiq_conf_values)
@pytest.mark.parametrize("only_ddis", only_ddis_values)
@pytest.mark.parametrize("remove_not_in_frame", remove_not_in_frame_values)
@pytest.mark.parametrize("only_divisible_by_three", only_divisible_by_three_values)
def test_nease_combinations(input, org, pred_ddi, p_value, delta, majiq_conf, only_ddis, remove_not_in_frame,
                  only_divisible_by_three):
    #print(f"input: {input}, org: {org}, pred_ddi: {pred_ddi}, p_value: {p_value}, delta: {delta}, ",
    #        f"majiq_conf: {majiq_conf}, only_ddis: {only_ddis}, remove_not_in_frame: {remove_not_in_frame}, ",
    #        f"only_divisible_by_three: {only_divisible_by_three}")
    database = input[0]
    table = readTable input
    events, info_tables, run_id = no.run_nease(table, org, {'db_type': database_type,
                                                                 'enrich_dbs': enrich_dbs,
                                                                 'p_value': p_value,
                                                                 'rm_not_in_frame': rm_not_in_frame,
                                                                 'divisible_by_3': divisible_by_3,
                                                                 'min_delta': min_delta,
                                                                 'majiq_confidence': majiq_confidence,
                                                                 'only_ddis': only_ddis,
                                                                 'confidences': confidences},
                                               input_data['splicing-events-file'].name,
                                               custom_name)