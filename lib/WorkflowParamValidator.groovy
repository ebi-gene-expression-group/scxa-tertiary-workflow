import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern TOKEN = Pattern.compile(/[A-Za-z0-9][A-Za-z0-9._+-]*/)
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern DECIMAL = Pattern.compile(/[0-9]+(\.[0-9]+)?/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Set TECHNOLOGIES = ['plate', 'droplet'] as Set
    private static final Set PUBLISH_MODES = ['copy', 'move', 'link', 'symlink', 'rellink'] as Set

    static void validate(def params) {
        requirePath(params, 'scanpy_scripts_container')
        requireEnum(params, 'technology', TECHNOLOGIES)
        optionalToken(params, 'batch_field')

        requirePath(params, 'dir_path')
        optionalPath(params, 'output_path')
        requirePath(params, 'result_dir_path')
        requireEnum(params, 'publish_dir_mode', PUBLISH_MODES)

        requireToken(params, 'representation')
        requireToken(params, 'celltype_field')
        requireList(params, 'neighbor_values', INTEGER)
        requireList(params, 'perplexity_values', DECIMAL)
        requireList(params, 'resolution_values', DECIMAL)
        requireToken(params, 'slotname')
        requireList(params, 'clustering_slotname', TOKEN)
        requireList(params, 'merged_group_slotname', TOKEN)
    }

    private static void requirePath(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", PATH_VALUE)
    }

    private static void optionalPath(def params, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "params.${name}", PATH_VALUE)
        }
    }

    private static void requireToken(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", TOKEN)
    }

    private static void optionalToken(def params, String name) {
        if (has(params, name) && params.get(name) != null && params.get(name).toString() != '') {
            assertPattern(params.get(name), "params.${name}", TOKEN)
        }
    }

    private static void requireEnum(def params, String name, Set allowed) {
        requireValue(params, "params.${name}", name)
        def text = params.get(name).toString()
        if (!(text in allowed)) {
            throw new IllegalArgumentException("params.${name} must be one of ${allowed}; got '${text}'")
        }
    }

    private static void requireList(def params, String name, Pattern pattern) {
        requireValue(params, "params.${name}", name)
        def values = params.get(name)
        def items = values instanceof Collection ? values : values.toString().split(',') as List
        if (items.isEmpty()) {
            throw new IllegalArgumentException("params.${name} must contain at least one value")
        }
        items.each { item ->
            assertPattern(item, "params.${name}", pattern)
        }
    }

    private static void requireValue(def params, String label, String key) {
        if (params == null || !has(params, key) || params.get(key) == null || params.get(key).toString() == '') {
            throw new IllegalArgumentException("Missing required workflow parameter ${label}")
        }
    }

    private static void assertPattern(value, String label, Pattern pattern) {
        def text = value == null ? '' : value.toString()
        if (!pattern.matcher(text).matches()) {
            throw new IllegalArgumentException("Invalid workflow parameter ${label}: '${text}'")
        }
    }

    private static boolean has(def params, String key) {
        params != null && params.containsKey(key)
    }
}
