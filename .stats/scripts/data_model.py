from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional


@dataclass
class ToolInfo:
    tool_id: str
    tool_name: str
    directory: str
    xml_files: list[str] = field(default_factory=list)
    categories: list[str] = field(default_factory=list)


@dataclass
class ToolEvent:
    tool_id: str
    tool_name: str
    directory: str
    event_type: str
    version: str
    old_version: Optional[str] = None
    commit_hash: str = ''
    commit_date: Optional[datetime] = None
    change_type: Optional[str] = None
    commit_message: str = ''

    def __post_init__(self):
        if isinstance(self.commit_date, str):
            self.commit_date = datetime.fromisoformat(self.commit_date)
