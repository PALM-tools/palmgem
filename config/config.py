import os
import sys
import yaml
import datetime
from collections import defaultdict

error_output = sys.stderr.write

def join(loader, node):
    seq = loader.construct_sequence(node)
    return ''.join([str(i) for i in seq])

## register the tag handler
yaml.add_constructor('!join', join)

def die(s, *args, **kwargs):
    """Write message to error output and exit with status 1."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')
    sys.exit(1)

def warn(s, *args, **kwargs):
    """Write message to error output."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')


class ConfigError(Exception):
    def __init__(self, desc, section=None, key=None):
        self.desc = desc
        self.section = section
        self.key = key

        # Build message
        s = ['Configuration error: ', desc]
        if section:
            s.extend([', item: ', ':'.join(section._get_path()+[key])])
            try:
                v = section._settings[key]
            except KeyError:
                s.append(', missing value')
            else:
                s.extend([', value=', str(v)])
        s.append('.')
        self.msg = ''.join(s)

    def __str__(self):
        return self.msg

class ConfigObj(object):
    """A recursive object within a hierarchical configuration, representing
    a (sub)section as a dictionary from the YAML configuration file. Child
    nodes may be accessed both by the dot notation (section.setting) and the
    item notation (section['setting']).
    """
    # We use __slots__ because we intend to hardly limit (and control) instance
    # members, so that we do not break many potential names of actual settings
    # that are accessed using the dot notation. For the same reason, member and
    # method names (mostly used internally anyway) start with an underscore.
    __slots__ = ['_parent', '_name', '_settings']

    def __init__(self, parent=None, name=None):
        self._parent = parent
        self._name = name
        self._settings = {}

    def __getattr__(self, name):
        try:
            return self._settings[name]
        except KeyError:
            raise AttributeError('Attribute {} not found. Possibly a missing '
                    'configuration setting in section {}.'.format(name,
                        ':'.join(self._get_path())))

    def __getitem__(self, key):
        try:
            return self._settings[key]
        except KeyError:
            raise KeyError('Key {} not found. Possibly a missing configuration '
                    'setting in section {}.'.format(key,
                        ':'.join(self._get_path())))

    def __iter__(self):
        return iter(self._settings.items())

    def _ingest_dict(self, d, overwrite=True, extend=True, check_exist=False):
        for k, v in d.items():
            if isinstance(v, ConfigObj):
                # we are actually ingesting a subtree - replace by its dict
                v = v._settings

            if isinstance(v, dict):
                # For a dictionary (top-level or with only dictionaries above,
                # i.e. a subsection), we recurse
                try:
                    vl = self._settings[k]
                except KeyError:
                    # not yet present: create a new empty child node
                    vl = ConfigObj(self, k)
                    self._settings[k] = vl
                try:
                    vl._ingest_dict(v, overwrite, extend, check_exist)
                except AttributeError:
                    raise ConfigError('Trying to replace a non-dictionary '
                            'setting with a dictionary', self, k)
            elif extend and isinstance(v, list):
                # We extend lists if requested
                vl = self._settings.setdefault(k, [])
                try:
                    vl.extend(v)
                except AttributeError:
                    raise ConfigError('Trying to extend a non-list setting with '
                            'a list', self, k)
            elif v is None and isinstance(self._settings.get(k), ConfigObj):
                # This is a special case: we are replacing an existing section
                # with None. That most probably means that the user has
                # presented an empty section (possibly with all values
                # commented out). In that case, we do not want to erase that
                # section. To actually erase a whole section, the user can
                # still present empty dictionary using the following syntax:
                # section_name: {}
                pass
            else:
                # Non-extended lists and all other objects are considered as
                # values and they are copied as-is (including their subtrees if
                # present). Non-null values are overwritten only if
                # overwrite=True.
                if overwrite:
                    if check_exist and k not in self._settings:
                        warn('WARNING: ignoring an unknown setting {}={}.',
                                ':'.join(self._get_path()+[k]), v)
                    self._settings[k] = v
                else:
                    if self._settings.get(k, None) is None:
                        self._settings[k] = v

    def _get_path(self):
        if self._parent is None:
            return []
        path = self._parent._get_path()
        path.append(self._name)
        return path


duration_units = {
        'd': 'days',
        'h': 'hours',
        'm': 'minutes',
        's': 'seconds',
        }

def parse_duration(section, item):
    def err():
        raise ConfigError('Bad specification of duration. The correct format is '
                '{num} {unit} [{num} {unit} ...], where {unit} is one of d, h, '
                'm, s. Example: "1 m 3.2 s".', section, item)

    try:
        s = section[item]
    except KeyError:
        err()

    words = s.split()
    n = len(words)
    if n % 2:
        err()

    d = defaultdict(int)
    for i in range(0, n, 2):
        ns, unit = words[i:i+2]
        try:
            num = int(ns)
        except ValueError:
            try:
                num = float(ns)
            except ValueError:
                err()
        try:
            u = duration_units[unit]
        except KeyError:
            err()
        d[u] += num

    return datetime.timedelta(**d)


def load_config(argv, cfg_default_path='config/default_config.yaml', cfg_default_path_share='config/default_share.yaml'):
    """Loads all configuration.

    Configuration is loaded in this order:
    1) initial configuration values
    2) configfile
    3) command-line options
    Each step may overwrite values from previous steps.
    """
    global cfg

    # load default configuration
    with open(cfg_default_path_share, 'r') as f:
        cfg._ingest_dict(yaml.load(f))

    # load default configuration
    with open(cfg_default_path, 'r') as f:
        cfg._ingest_dict(yaml.load(f))

    # load settings from user configfile (if available)
    if argv.config:
        with open('config/' + argv.config, 'r') as config_file:
            cfg._ingest_dict(yaml.load(config_file), check_exist=True)

    # extras based on the specific case
    case_schema = cfg.domain.name if cfg.domain.scenario == "" else cfg.domain.name + '_' + cfg.domain.scenario
    cfg.domain._settings['case_schema'] = case_schema

    # static driver netcdf file name
    if not 'static_driver_file' in cfg.domain._settings.keys():
        if cfg.domain.scenario == "":
            static_driver_file = cfg.domain.name + '_static.nc'
        else:
            static_driver_file = cfg.domain.name + '_' + cfg.domain.scenario + '_static.nc'
        cfg.domain._settings['static_driver_file'] = os.path.join('output', static_driver_file)

    # name file of visual check according to case_schema
    cfg.visual_check._settings['path'] = os.path.join('visual_check', case_schema)

    if not os.path.isdir(os.path.join('visual_check')):
        os.mkdir(os.path.join('visual_check'))

    # create file for visual check
    if cfg.visual_check.enabled and not os.path.isdir(cfg.visual_check.path):
        os.mkdir(cfg.visual_check.path)

    #create file for netcdf outputs
    if not os.path.isdir('output'):
        os.mkdir('output')

    # name log file according to case_schema
    if not os.path.isdir('logs'):
        os.mkdir('logs')

    cfg.logs._settings['path'] = os.path.join('logs', case_schema)

    # check log level, if sublevel are not specified, set them as general log level
    for lvl in cfg.logs._settings.keys():
        if 'level_' in lvl:
            if cfg.logs[lvl] == -1:
                cfg.logs._settings[lvl] = cfg.logs.level


cfg = ConfigObj()