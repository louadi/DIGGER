// Get all elements with the class "span_value"
const spanElements = document.querySelectorAll('.example_value');

// Loop through each span element and add event listener
spanElements.forEach(function(spanElement) {
    spanElement.addEventListener('click', function() {
        // Get the value of the data-input-field attribute
        const inputFieldId = spanElement.getAttribute('data-input-field');
        // Get the corresponding input field
        const inputField = document.getElementById(inputFieldId);
        // Set the value of the input field to the text content of the clicked span element
        inputField.value = spanElement.textContent;
    });
});