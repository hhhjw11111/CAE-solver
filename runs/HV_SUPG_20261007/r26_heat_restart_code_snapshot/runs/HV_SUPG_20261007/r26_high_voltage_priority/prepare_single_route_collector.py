"""Prepare one collector for the one retained high-voltage continuation."""
from pathlib import Path
import ast

here=Path(__file__).absolute().parent
text=(here/'collect_high_voltage.py').read_text()
old="SOURCES={'full_heat_voltage':BASE/'r25_high_voltage_followthrough/target_V080',\n         'trust_heat_8V':BASE/'r26_log_A64_trust_heat_V080'}"
new="SOURCES={'trust_heat_restart_8V':BASE/'r26_log_A64_trust_heat_restart_V080'}"
assert old in text
text=text.replace(old,new).replace("OUT=HERE/'results'","OUT=HERE/'results_restart'")
marker="    if (existing/'execution.json').exists():add_record(6.25,existing,'full_heat_voltage')"
assert marker in text
text=text.replace(marker,marker+"\n    add_record(6.375,HERE/'results/full_heat_voltage_V0638','completed_fixed_voltage_replay')")
ast.parse(text)
dest=here/'collect_high_voltage_restart.py'
with dest.open('x') as f:f.write(text)
print(dest)
