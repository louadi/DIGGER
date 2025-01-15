import pytest
import pandas as pd
# Set up according to https://djangostars.com/blog/django-pytest-testing/ and
# https://stackoverflow.com/questions/60369047/pytest-using-parametized-fixture-vs-pytest-mark-parametrize

from .Process import nease_output as no
# Run with the command pytest, add -s if you want prints (might not work now, as we added pytest-xdist), add -n auto if you want to run in parallel

location="./domain"
input_files = [("Standard", "human", location+"/test_files/Standard.tsv"),
               ("MAJIQ", "human", location+"/test_files/MAJIQ.tsv"),
               ("Whippet", "mouse", location+"/test_files/Whippet.diff"),
               ("rmats", "mouse", location+"/test_files/rMATS.txt")]

predicted_ddis = [[], ["high"], ["mid"], ["low"], ["high", "mid"], ["mid", "low"],  ["high", "low"], ["high", "mid", "low"]]
p_values = [1, 0.05, 0]
deltas = [1, 0.05, -1]
majiq_conf_values = [1, 0.95, 0]
only_ddis_values = [True, False]
remove_not_in_frame_values = [True, False]
only_divisible_by_three_values = [True, False]


# Use pytest.mark.parametrize to generate all combinations of the parameters
@pytest.mark.parametrize("input", input_files)
@pytest.mark.parametrize("pred_ddi", predicted_ddis)
@pytest.mark.parametrize("p_value", p_values)
@pytest.mark.parametrize("delta", deltas)
@pytest.mark.parametrize("majiq_conf", majiq_conf_values)
@pytest.mark.parametrize("only_ddis", only_ddis_values)
@pytest.mark.parametrize("remove_not_in_frame", remove_not_in_frame_values)
@pytest.mark.parametrize("only_divisible_by_three", only_divisible_by_three_values)
@pytest.mark.django_db
def test_nease_combinations(input, pred_ddi, p_value, delta, majiq_conf, only_ddis, remove_not_in_frame,
                  only_divisible_by_three):
    #print(f"input: {input}, org: {org}, pred_ddi: {pred_ddi}, p_value: {p_value}, delta: {delta}, ",
    #        f"majiq_conf: {majiq_conf}, only_ddis: {only_ddis}, remove_not_in_frame: {remove_not_in_frame}, ",
    #        f"only_divisible_by_three: {only_divisible_by_three}")
    database = input[0]
    org = input[1]
    # If the database is not MAJIQ, dont run everything three times but only for one confidence value
    if database != "MAJIQ" and majiq_conf != 0.95:
        return
    filepath = input[2]
    if database == "Whippet":
        table = no.read_extra_spaces(filepath)
    else:
        table = pd.read_table(filepath)
    try:
        events, info_tables, run_id = no.run_nease(table, org, {'db_type': database,
                                                                'p_value': p_value,
                                                                'rm_not_in_frame': remove_not_in_frame,
                                                                'divisible_by_3': only_divisible_by_three,
                                                                'min_delta': delta,
                                                                'majiq_confidence': majiq_conf,
                                                                'only_ddis': only_ddis,
                                                                'confidences': pred_ddi},
                                                   "filename",
                                                   "custom_name")
    except ValueError as e:
        msg = e.args[0]
        if "No significant events" in msg:
            return
        else:
            raise e