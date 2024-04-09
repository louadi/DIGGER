// Derived from here https://stackoverflow.com/questions/45588138/typeahead-js-in-django-project

var gene_symbol = new Bloodhound({
            datumTokenizer: Bloodhound.tokenizers.obj.whitespace('gene_symbol'),
            queryTokenizer: Bloodhound.tokenizers.whitespace,
            prefetch: document.body.getAttribute('data-typeahead') + '?q=%QUERY&o=%ORGANISM',
            remote: {
                url: document.body.getAttribute('data-typeahead') + '?q=%QUERY&o=%ORGANISM',
                wildcard: '%QUERY',
                replace: function(url, query) {
                    var organism = $('#org_select').val();
                    // only consider the last characters after any comma if there is one (support for multi input)
                    var last_query = query.split(',').pop().trim();
                    if (last_query.length > 0) {
                        query = last_query;
                    }
                    return url.replace('%QUERY', query).replace('%ORGANISM', organism);
                }
            }
 });

$('#remote .typeahead').typeahead(null, {
    name: 'gene-symbol',
    display: 'gene_symbol',
    source: gene_symbol,
    templates: {
        // Derived from here https://github.com/twitter/typeahead.js/issues/1031#issuecomment-66567928
        suggestion: function (data) {
            return '<p><strong>' + data.gene_symbol + '</strong> - ' + data.ensembl_id + '</p>';
        }
    }

});